#include <algorithm>

#include <FieldBasis.H>
#include <EXPException.H>
#include <DiskModels.H>
#include <YamlCheck.H>
#include <exputils.H>
#include <gaussQ.H>

#ifdef HAVE_FE_ENABLE
#include <cfenv>
#endif

extern double plgndr(int l, int m, double x);

static double Ylm(int l, int m)
{
  m = abs(m);
  return sqrt( (2.0*l+1)/(4.0*M_PI) ) *
    exp(0.5*(lgamma(1.0+l-m) - lgamma(1.0+l+m)));
}


namespace BasisClasses
{
  const std::set<std::string>
  FieldBasis::valid_keys = {
    "modelname",
    "dof",
    "rmapping",
    "rmin",
    "rmax",
    "ascl",
    "delta",
    "lmax",
    "nmax",
    "model"
  };

  void FieldBasis::addPSFunction(FieldBasis::PSFunction func,
				 std::vector<std::string>& labels)
  {
    // Test return dimensionality
    //
    double z = 0.01;
    PS3 pos{z, z, z}, vel{z, z, z};
    auto p = func(z, pos, vel);
    if (p.size() != labels.size()) {
      std::ostringstream sout;
      sout << "FieldBasis::register mismatch between field dimension <"
	   << p.size() << "> and label dimension <" << labels.size() << ">";
      throw std::runtime_error(sout.str());
    }

    // Allocate coefficient storage
    //
    nfld = p.size() + 2;
    allocateStore();

    // Okay to register
    //
    fieldFunc = func;
    for (auto & v : labels) fieldLabels.push_back(v);
  }

  // Deploy instance based on configuration file 
  //
  void FieldBasis::configure()
  {
    nfld      = 2;		// Weight and density fields by default
    lmax      = 4;
    nmax      = 10;
    rmin      = 1.0e-4;
    rmax      = 2.0;
    ascl      = 0.01;
    delta     = 0.005;
    rmapping  = 0.05;
    dof       = 3;
    model     = "file";
    name      = "field";
    modelname = "SLGridSph.model";

    initialize();

    // Check for valid model type
    //
    std::vector<std::string> modelTypes { "file", "expon" };
    auto it = std::find(modelTypes.begin(), modelTypes.end(), model);
    if (it == modelTypes.end()) {
      std::ostringstream sout;
      sout << "FieldBasis: found type <" << model << ">.  Must be one of ";
      for (auto v : modelTypes) sout << " " << v;
      throw std::runtime_error(sout.str());
    }
    
    // Check dof value
    //
    if (dof!=2 and dof!=3) {
      std::ostringstream sout;
      sout << "FieldBasis: found " << dof << " for dof.  Must be 2 or 3.";
      throw std::runtime_error(sout.str());
    }

    // Assign coordinate type
    //
    if (dof==2) coordinates = Coord::Cylindrical;
    else        coordinates = Coord::Spherical;

    // Allocate storage for coefficients
    //
    nt = omp_get_max_threads();
    store.resize(nt);
    coefs.resize(nt);
    massT.resize(nt);
    usedT.resize(nt);
    
    // Create model needed for density prefactor in OrthoFunction
    //
    if (model == "file") {
      std::vector<double> r, d;
      std::ifstream in(modelname);
      if (not in) throw std::runtime_error("Error opening file: " + modelname);
    
      std::string line;
      while (std::getline(in, line)) {
	auto pos = line.find_first_of("!#");
	if (pos == std::string::npos) {
	  std::istringstream iss(line);
	  double x, y;
	  iss >> x >> y;
	  if (iss) {
	    r.push_back(x);
	    d.push_back(y);
	  }
	}
      }
      
      // Compute interpolation functionoid
      //
      double rmin = r.front(), rmax = r.back();
      
      interp = std::make_shared<Linear1d>(r, d);
      densfunc = [this](double r)
      {
	return this->interp->eval(r);
      };
      
    } else if (model == "expon") {
      
      // Density functionoid for the exponential
      //
      densfunc = [&](double r)
      {
	return exp(-r/ascl) *
	  0.5*(1.0 + std::erf((rmax - 5.0*delta - r)/delta)) / ascl;
      };
    
    } else {
      throw InternalError("FieldBasis: model logic failure?! "
			  "You should not be here...",
			  __FILE__, __LINE__);
    }
    
    // Generate the orthogonal function instance
    //
    ortho = std::make_shared<OrthoFunction>(nmax, densfunc, rmin, rmax, rmapping, dof);

    // Initialize fieldlabels
    //
    fieldLabels.clear();
    fieldLabels.push_back("weight");
    fieldLabels.push_back("density");

    // Debug
    //
    if (true) {
      auto tst = ortho->testOrtho();
      double worst = 0.0;
      for (int i=0; i<tst.rows(); i++) {
	for (int j=0; j<tst.rows(); j++) {
	  if (i==j) worst = std::max<double>(worst, fabs(1.0 - tst(i, j)));
	  else      worst = std::max<double>(worst, fabs(tst(i, j)));
	}
      }
      if (myid==0) {
	std::cout << "FieldBasis::orthoTest: worst=" << worst << std::endl;
	ortho->dumpOrtho("fieldbasis_ortho.dump");
      }
    }
  }

  void FieldBasis::allocateStore() 
  {
    for (int t=0; t<nt; t++) {
      if (dof==2) {
	store[t].resize(nfld*(lmax+1)*nmax);
	coefs[t] = std::make_shared<coefType>(store[t].data(), nfld, lmax+1, nmax);
      } else {
	store[t].resize(nfld*(lmax+1)*(lmax+2)/2*nmax);
	coefs[t] = std::make_shared<coefType>(store[t].data(), nfld, (lmax+1)*(lmax+2)/2, nmax);
      }
    }
  }

  void FieldBasis::initialize()
  {
    // Remove matched keys
    //
    // for (auto v : valid_keys) current_keys.erase(v);
  
    // Assign values from YAML
    //
    try {
      if (conf["modelname"])	modelname = conf["modelname"].as<std::string>();
      if (conf["model"    ])    model     = conf["model"    ].as<std::string>();
      if (conf["nfld"     ])    nfld      = conf["nfld"     ].as<int>();
      if (conf["lmax"     ])    lmax      = conf["lmax"     ].as<int>();
      if (conf["nmax"     ])    nmax      = conf["nmax"     ].as<int>();
      if (conf["dof"      ])    dof       = conf["dof"      ].as<int>();
      if (conf["rmin"     ])    rmin      = conf["rmin"     ].as<double>();
      if (conf["rmax"     ])    rmax      = conf["rmax"     ].as<double>();
      if (conf["ascl"     ])    ascl      = conf["ascl"     ].as<double>();
      if (conf["delta"    ])    delta     = conf["delta"    ].as<double>();
      if (conf["rmapping" ])    rmapping  = conf["rmapping" ].as<double>();
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameters in FieldBasis: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      throw std::runtime_error("FieldBasis::initialize: error parsing YAML");
    }
  }

  void FieldBasis::load_coefs(CoefClasses::CoefStrPtr cf, double time)
  {
    // Assign the coeffiicient container and the time
    // 
    coefret = cf;
    coefret->time = time;

    std::ostringstream sout; sout << node;

    // Copy the data structure
    if (dof==2) {
      auto p = dynamic_pointer_cast<CoefClasses::CylFldStruct>(cf);
      p->nfld  = nfld;
      p->mmax  = lmax;
      p->nmax  = nmax;
      p->allocate();
      p->store = store[0];
      p->buf   = sout.str();
    } else {
      auto p = dynamic_pointer_cast<CoefClasses::SphFldStruct>(cf);
      p->nfld  = nfld;
      p->lmax  = lmax;
      p->nmax  = nmax;
      p->allocate();
      p->store = store[0];
      p->buf   = sout.str();
    }
  }

  void FieldBasis::reset_coefs()
  {
    used = 0;
    totalMass = 0.0;
    for (auto & v : usedT) v = 0;
    for (auto & v : massT) v = 0.0;
    for (auto & v : store) v.setZero();
  }

  void FieldBasis::set_coefs(CoefClasses::CoefStrPtr c)
  {
    if (dof==2) {
      auto p = dynamic_cast<CoefClasses::CylFldStruct*>(c.get());
      coefs[0] = p->coefs;
      store[0] = p->store;

      // Sanity test dimensions
      if (nfld!=p->nfld || lmax!=p->mmax || nmax!=p->nmax) {
	std::ostringstream serr;
	serr << "FieldBasis::set_coefs: dimension error! "
	     << " nfld [" << nfld << "!= " << p->nfld << "]"
	     << " mmax [" << lmax << "!= " << p->mmax << "]"
	     << " nmax [" << nmax << "!= " << p->nmax << "]";
	throw std::runtime_error(serr.str());
      }

    } else {
      auto p = dynamic_cast<CoefClasses::SphFldStruct*>(c.get());
      coefs[0] = p->coefs;
      store[0] = p->store;

      // Sanity test dimensions
      if (nfld!=p->nfld || lmax!=p->lmax || nmax!=p->nmax) {
	std::ostringstream serr;
	serr << "FieldBasis::set_coefs: dimension error! "
	     << " nfld [" << nfld << "!= " << p->nfld << "]"
	     << " lmax [" << lmax << "!= " << p->lmax << "]"
	     << " nmax [" << nmax << "!= " << p->nmax << "]";
	throw std::runtime_error(serr.str());
      }

    }
  }

  void FieldBasis::accumulate(double mass,
			      double x, double y, double z,
			      double u, double v, double w)
  {
    constexpr std::complex<double> I(0, 1);
    constexpr double fac0 = 1.0/sqrt(4*M_PI);

    int tid = omp_get_thread_num();
    PS3 pos{x, y, z}, vel{u, v, w};

    // Compute the field value array
    //
    std::vector<double> vec;
    if (fieldFunc) vec = fieldFunc(mass, pos, vel);
    
    // Compute spherical/polar coordinates
    //
    double R   = sqrt(x*x + y*y);
    double r   = sqrt(R*R + z*z);
    double phi = atan2(y, x);
    double cth = z/(r + 1.0e-18);
    
    usedT[tid] ++;		// Count used particles
    massT[tid] += mass;		// Sum particles' mass

    if (dof==2) {
      
      auto p = (*ortho)(R);
    
      (*coefs[tid])(0, 0, 0) += mass*p(0)*fac0;

      for (int m=0; m<=lmax; m++) {
	
	std::complex<double> P = std::exp(-I*(phi*m));
	
	for (int n=0; n<nmax; n++) {

	  (*coefs[tid])(1, m, n) += mass*P*p(n);

	  for (int k=0; k<vec.size(); k++)
	    (*coefs[tid])(k+2, m, n) += mass*P*p(n)*vec[k];
	}
      }	 
      
    } else {
    
      auto p = (*ortho)(r);
      
      (*coefs[tid])(0, 0, 0) += mass*p(0);

      for (int l=0, lm=0; l<=lmax; l++) {

	double s = 1.0;

	for (int m=0; m<=l; m++, lm++) {
	  
	  // Spherical harmonic value
	  std::complex<double> P =
	    std::exp(-I*(phi*m))*Ylm(l, m)*plgndr(l, m, cth) * s;

	  s *= -1.0;		// Flip sign for next evaluation

	  for (int n=0; n<nmax; n++) {

	    (*coefs[tid])(1, lm, n) += mass*P*p(n);

	    for (int k=0; k<vec.size(); k++)
	      (*coefs[tid])(k+2, lm, n) += mass*P*p(n)*vec[k];
	  }
	}
      }
    }
  }

  void FieldBasis::make_coefs()
  {
    // Sum over threads
    //
    for (int t=1; t<coefs.size(); t++) {
      store[0] += store[t];
      usedT[0] += usedT[t];
      massT[0] += massT[t];
    }
      
    if (use_mpi) {
      MPI_Allreduce(MPI_IN_PLACE, store[0].data(), store[0].size(),
		    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(&usedT[0], &used,      1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&massT[0], &totalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } else {
      used = usedT[0];
      totalMass = massT[0];
    }

    // Create a structure of the appropriate type, if needed
    //
    if (not coefret) {
      if (dof==2) coefret = std::make_shared<CoefClasses::CylFldStruct>();
      else        coefret = std::make_shared<CoefClasses::SphFldStruct  >();
				// Add the configuration data
      std::ostringstream sout; sout << node;
      coefret->buf = sout.str();
    }

    // Assign the data by downcasting.  Ideally, the 'assign' method
    // would be in the base class but the signature depends on data
    // type in the derived class.
    //
    if (dof==2) {
      auto p = dynamic_pointer_cast<CoefClasses::CylFldStruct>(coefret);
      p->assign(store[0], nfld, lmax, nmax);
    } else {
      auto p = dynamic_pointer_cast<CoefClasses::SphFldStruct>(coefret);
      p->assign(store[0], nfld, lmax, nmax);
    }
  }

  //! Retrieve the coefficients 
  CoefClasses::CoefStrPtr FieldBasis::getCoefficients()
  {
    if (dof==2) {
      auto cstr = std::make_shared<CoefClasses::CylFldStruct>();
      cstr->assign(store[0], nfld, lmax, nmax);
      return cstr;
    } else {
      auto cstr = std::make_shared<CoefClasses::SphFldStruct>();
      cstr->assign(store[0], nfld, lmax, nmax);
      return cstr;
    }
  }
  
  std::vector<double>
  FieldBasis::sph_eval(double r, double costh, double phi)
  {
    constexpr std::complex<double> I(0, 1);

    if (dof==2) {
      std::vector<double> ret(nfld, 0);
      auto p = (*ortho)(r*sqrt(fabs(1.0 - costh*costh)));
      for (int i=0; i<nfld; i++) {
	for (int m=0; m<=lmax; m++) {
	  double fac = m>0 ? 2.0 : 1.0;
	  for (int n=0; n<nmax; n++) {
	    ret[i] += fac *
	      std::real((*coefs[0])(i, m, n)*exp(I*(phi*m)))*p[n];
	  }
	}
      }
      return ret;
    } else {
      static bool firstime = true;
      if (firstime) {
	int tid = omp_get_thread_num();
	std::ostringstream file;
	file << "field.bin." << tid;
	std::ofstream fout(file.str());
	const auto& d = coefs[0]->dimensions();
	fout << "Dim size: " << d.size();
	for (int i=0; i<d.size(); i++) fout << ", dim " << i << ": " << d[i];
	fout << std::endl;
	for (auto v : store[0]) fout << v << std::endl;
	firstime = false;
      }

      std::vector<double> ret(nfld, 0);
      auto p = (*ortho)(r);
      for (int i=0; i<nfld; i++) {
	for (int l=0, L=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++, L++) {
	    double fac = m>0 ? 2.0 : 1.0;
	    for (int n=0; n<nmax; n++) {
	      ret[i] += std::real((*coefs[0])(i, L, n)*exp(I*(phi*m)))*p[n]*
		Ylm(l, m)*plgndr(l, m, costh) * fac;
	    }
	  }
	}
      }
      return ret;
    }
  }

  std::vector<double>
  FieldBasis::cyl_eval(double R, double z, double phi)
  {
    double r = sqrt(R*R + z*z);
    double costh = z/(r + 1.0e-18);

    return sph_eval(r, costh, phi);
  }

  std::vector<double>
  FieldBasis::crt_eval(double x, double y, double z)
  {
    double r = sqrt(x*x + y*y + z*z);
    double costh = z/(r + 1.0e-18), phi = atan2(y, x);
    
    return sph_eval(r, costh, phi);
  }

  FieldBasis::BasisArray
  FieldBasis::getBasis(double logxmin, double logxmax, int numgrid)
  {
    // Assign return storage
    //
    BasisArray ret(nmax, numgrid);

    // Radial grid spacing
    double dx = (logxmax - logxmin)/numgrid;

    for (int i=0; i<numgrid; i++) {
      double R = pow(10.0, logxmin + dx*i);
      ret.col(i) = (*ortho)(R);
    }
    
    return ret;
  }

  // Generate coeffients from a particle reader
  CoefClasses::CoefStrPtr FieldBasis::createFromReader
  (PR::PRptr reader, std::vector<double> ctr)
  {
    CoefClasses::CoefStrPtr coef;
    
    if (dof==3)
      coef = std::make_shared<CoefClasses::SphFldStruct>();
    else if (dof==2)
      coef = std::make_shared<CoefClasses::CylFldStruct>();
    else {
      std::ostringstream sout;
      sout << "FieldBasis::createCoefficients: dof must be 2 or 3"
	   << std::endl;
      throw std::runtime_error(sout.str());
    }
    
    // Is center non-zero?
    //
    bool addCenter = false;
    for (auto v : ctr) {
      if (v != 0.0) addCenter = true;
    }

    // Add the expansion center metadata
    //
    if (addCenter) coef->ctr = ctr;

    std::vector<double> pp(3), vv(3);

    reset_coefs();
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {

      bool use = false;
      
      if (ftor) {
	pp.assign(p->pos, p->pos+3);
	vv.assign(p->vel, p->vel+3);
	use = ftor(p->mass, pp, vv, p->indx);
      } else {
	use = true;
      }

      if (use) accumulate(p->mass,
			  p->pos[0]-ctr[0],
			  p->pos[1]-ctr[1],
			  p->pos[2]-ctr[2],
			  p->vel[0],
			  p->vel[1],
			  p->vel[2]);

    }
    make_coefs();
    load_coefs(coef, reader->CurrentTime());
    return coef;
  }

  // Generate coefficients from a phase-space table
  void FieldBasis::initFromArray(std::vector<double> ctr)
  {
    if (dof==3)
      coefret = std::make_shared<CoefClasses::SphFldStruct>();
    else if (dof==2)
      coefret = std::make_shared<CoefClasses::CylFldStruct>();
    else {
      std::ostringstream sout;
      sout << "FieldBasis::createCoefficients: dof must be 2 or 3"
	   << std::endl;
      throw std::runtime_error(sout.str());
    }
      
    // Is center non-zero?
    //
    bool addCenter = false;
    for (auto v : ctr) {
      if (v != 0.0) addCenter = true;
    }

    // Add the expansion center metadata
    //
    if (addCenter) coefret->ctr = ctr;
    
    // Register the center
    //
    coefctr = ctr;

    // Clean up for accumulation
    //
    reset_coefs();
    coefindx = 0;
  }

  // Accumulate coefficient contributions from arrays
  void FieldBasis::addFromArray(Eigen::VectorXd& m, RowMatrixXd& p,
				bool roundrobin, bool posvelrows)
  {
    // Sanity check: is coefficient instance created?  This is not
    // foolproof.  It is really up the user to make sure that a call
    // to initFromArray() comes first.
    //
    if (not coefret) {
      std::string msg =
	"FieldBasis::addFromArray: you must initialize coefficient accumulation "
	"with a call to FieldBasis::initFromArray()";
      throw std::runtime_error(msg);
    }

    std::vector<double> p1(3), v1(3);

    if (posvelrows) {

      if (p.rows()<6) {

	std::ostringstream msg;
	msg << "Basis::addFromArray: you must pass a position array with at "
	  "least three six for x, y, z, u, v, w.  Yours has " << p.rows() << ".";
	throw std::runtime_error(msg.str());
      }

      for (int n=0; n<p.cols(); n++) {

	if (n % numprocs==myid or not roundrobin) {

	  bool use = true;
	  if (ftor) {
	    for (int k=0; k<3; k++) {
	      p1[k] = p(k+0, n);
	      v1[k] = p(k+3, n);
	    }
	    use = ftor(m(n), p1, v1, coefindx);
	  } else {
	    use = true;
	  }
	  coefindx++;
	  
	  if (use) accumulate(m(n),
			      p(0, n)-coefctr[0],
			      p(1, n)-coefctr[1],
			      p(2, n)-coefctr[2],
			      p(3, n),
			      p(4, n),
			      p(5, n));
	}
      }
      
    } else {

      if (p.cols()<6) {
	std::ostringstream msg;
	msg << "Basis::addFromArray: you must pass a position array with at "
	  "least six columns for x, y, z, u, v, w.  Yours has " << p.cols() << ".";
	throw std::runtime_error(msg.str());
      }

      for (int n=0; n<p.rows(); n++) {

	if (n % numprocs==myid or not roundrobin) {

	  bool use = true;
	  if (ftor) {
	    for (int k=0; k<3; k++) {
	      p1[k] = p(n, k+0);
	      v1[k] = p(n, k+3);
	    }
	    use = ftor(m(n), p1, v1, coefindx);
	  } else {
	    use = true;
	  }
	  coefindx++;
	  
	  if (use) accumulate(m(n),
			      p(n, 0)-coefctr[0],
			      p(n, 1)-coefctr[1],
			      p(n, 2)-coefctr[2], 
			      p(n, 3),
			      p(n, 4),
			      p(n, 5));
	}
      }
    }
  }

  std::vector<double> cylVel(double mass,
			     FieldBasis::PS3& pos, FieldBasis::PS3& vel)
  {
    auto [x, y, z] = pos;
    auto [u, v, w] = vel;
    
    double R   = sqrt(x*x + y*y) + 1.0e-18;
    double vr  = (u*x + v*y)/R;
    double vp  = (u*y - v*x)/R;
    
    return {vr, w, vp, vr*vr, w*w, vp*vp};
  }

  std::vector<double> sphVel(double mass,
			     FieldBasis::PS3& pos, FieldBasis::PS3& vel)
  {
    auto [x, y, z] = pos;
    auto [u, v, w] = vel;
    
    double R   = sqrt(x*x + y*y) + 1.0e-18;
    double r   = sqrt(R*R + z*z);
    
    double vr  = (u*x + v*y + w*z)/r;
    double vt  = (u*z*x + v*z*y - w*R)/R/r;
    double vp  = (u*y - v*x)/R;

    return {vr, vt, vp, vr*vr, vt*vt, vp*vp};
  }

  std::vector<double> crtVel(double mass,
			     FieldBasis::PS3& pos, FieldBasis::PS3& vel)
  {
    auto [u, v, w] = vel;
    return {u, v, w, u*u, v*v, w*w};
  }


  void VelocityBasis::assignFunc()
  {
    fieldLabels.clear();

    // Field labels (force field labels added below)
    //
    fieldLabels.push_back("weight");
    fieldLabels.push_back("density");

    if (coordinates == Coord::Cylindrical) {
      fieldLabels.push_back("v_R");
      fieldLabels.push_back("v_z");
      fieldLabels.push_back("v_p");
      fieldLabels.push_back("v_R^2");
      fieldLabels.push_back("v_z^2");
      fieldLabels.push_back("v_p^2");
      fieldFunc = cylVel;
    } else if (coordinates == Coord::Cartesian) {
      fieldLabels.push_back("v_x");
      fieldLabels.push_back("v_y");
      fieldLabels.push_back("v_z");
      fieldLabels.push_back("v_x^2");
      fieldLabels.push_back("v_y^2");
      fieldLabels.push_back("v_z^2");
      fieldFunc = crtVel;
    } else if (coordinates == Coord::None) {
      fieldLabels.push_back("v_x");
      fieldLabels.push_back("v_y");
      fieldLabels.push_back("v_z");
      fieldLabels.push_back("v_x^2");
      fieldLabels.push_back("v_y^2");
      fieldLabels.push_back("v_z^2");
      fieldFunc = crtVel;
    } else {
      fieldLabels.push_back("v_r");
      fieldLabels.push_back("v_t");
      fieldLabels.push_back("v_p");
      fieldLabels.push_back("v_r^2");
      fieldLabels.push_back("v_t^2");
      fieldLabels.push_back("v_p^2");
      fieldFunc = sphVel;
    }

    // Allocate storage
    //
    nfld = fieldLabels.size();
    allocateStore();
  }

  VelocityBasis::VelocityBasis(const YAML::Node& conf) :
    FieldBasis(conf, "VelocityBasis")
  {
    name = "velocity";
    assignFunc();
  }

  VelocityBasis::VelocityBasis(const std::string& confstr) :
    FieldBasis(confstr, "VelocityBasis")
  {
    name = "velocity";
    assignFunc();
  }


}
// END namespace BasisClasses
