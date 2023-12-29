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

extern double Ylm01(int ll, int mm);
extern double plgndr(int l, int m, double x);

namespace BasisClasses
{
  const std::set<std::string>
  FieldBasis::valid_keys = {
    "filename",
    "name",
    "dof",
    "scale",
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

    // Okay to register
    //
    fieldFunc   = func;
    fieldLabels = labels;
  }

  // Deploy instance based on configuration file 
  //
  void FieldBasis::configure()
  {
    lmax     = 4;
    nmax     = 10;
    rmin     = 1.0e-4;
    rmax     = 2.0;
    ascl     = 0.01;
    delta    = 0.005;
    scale    = 0.05;
    dof      = 3;
    model    = "file";
    filename = "SLGridSph.model";

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
    int nt = omp_get_max_threads();
    store.resize(nt);
    coefs.resize(nt);
    massT.resize(nt);
    usedT.resize(nt);
    
    for (int t=0; t<nt; t++) {
      if (dof==2) {
	store[t].resize(3*(lmax+1)*nmax);
	coefs[t] = std::make_shared<coefType>(store[t].data(), 3, lmax+1, nmax);
      } else {
	store[t].resize(4*(lmax+1)*(lmax+2)/2*nmax);
	coefs[t] = std::make_shared<coefType>(store[t].data(), 4, (lmax+1)*(lmax+2)/2, nmax);
      }
    }
       
    // Create model needed for density prefactor in OrthoFunction
    //
    if (model == "file") {
      std::vector<double> r, d;
      std::ifstream in(filename);
      if (not in) throw std::runtime_error("Error opening file: " + filename);
    
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
    ortho = std::make_shared<OrthoFunction>(nmax, densfunc, rmin, rmax, scale);
  }

  void FieldBasis::initialize()
  {
    // Remove matched keys
    //
    // for (auto v : valid_keys) current_keys.erase(v);
  
    // Assign values from YAML
    //
    try {
      if (conf["filename"])	filename = conf["filename"].as<std::string>();
      if (conf["lmax"])         lmax     = conf["lmax" ].as<int>();
      if (conf["nmax"])         nmax     = conf["nmax" ].as<int>();
      if (conf["dof"])          dof      = conf["dof"  ].as<int>();
      if (conf["rmin"])         rmin     = conf["rmin" ].as<double>();
      if (conf["rmax"])         rmax     = conf["rmax" ].as<double>();
      if (conf["ascl"])         ascl     = conf["ascl" ].as<double>();
      if (conf["delta"])        delta    = conf["delta"].as<double>();
      if (conf["scale"])        scale    = conf["scale"].as<double>();
      if (conf["model"])        model    = conf["model"].as<std::string>();
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameters in FieldBasis: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      if (use_mpi) MPI_Finalize();
      exit(-1);
    }
  }

  void FieldBasis::load_coefs(CoefClasses::CoefStrPtr cf, double time)
  {
    coefstr = cf;		// Assign the coeffiicient container
    coefstr->time = time;	// and the time

    std::ostringstream sout; sout << node;

    // Copy the data structure
    if (dof==2) {
      auto p = dynamic_pointer_cast<CoefClasses::PolarVelStruct>(cf);
      p->mmax  = lmax;
      p->nmax  = nmax;
      p->allocate();
      p->store = store[0];
      p->buf   = sout.str();
    } else {
      auto p = dynamic_pointer_cast<CoefClasses::SphVelStruct>(cf);
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
      auto p = dynamic_cast<CoefClasses::PolarVelStruct*>(c.get());
      coefs[0] = p->coefs;
    } else {
      auto p = dynamic_cast<CoefClasses::SphVelStruct*>(c.get());
      coefs[0] = p->coefs;
    }
  }

  void FieldBasis::accumulate(double mass,
			      double x, double y, double z,
			      double u, double v, double w)
  {
    constexpr std::complex<double> I(0, 1);
    int tid = omp_get_thread_num();
    PS3 pos{x, y, z}, vel{u, v, w};

    auto vec = fieldFunc(mass, pos, vel);
    
    double R   = sqrt(x*x + y*y);
    double r   = sqrt(R*R + z*z);
    double phi = atan2(y, x);
    double cth = z/(r + 1.0e-18);
    
    if (dof==2) {
      
      auto p = (*ortho)(R);
    
      for (int m=0; m<=lmax; m++) {
	std::complex<double> P = std::exp(I*(phi*m));
	for (int n=0; n<nmax; n++) {
	  (*coefs[tid])(0, m, n) += mass*P*p(n);
	  for (int k=0; k<p.size(); k++)
	    (*coefs[tid])(k+1, m, n) += mass*P*p(n)*vec[k];
	}
      }	 
      
    } else {
    
      auto p = (*ortho)(r);
      
      for (int l=0, k=0; l<=lmax; l++) {
	for (int m=0; m<=l; m++, k++) {
	  std::complex<double> P =
	    std::exp(I*(phi*m))*Ylm01(l, m)*plgndr(l, m, cth);
	
	  for (int n=0; n<nmax; n++) {
	    (*coefs[tid])(0, m, n) += mass*P*p(n);
	    for (int k=0; k<p.size(); k++)
	      (*coefs[tid])(k+1, m, n) += mass*P*p(n)*vec[k];
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
      *coefs[0] += *coefs[t];
      usedT[0] += usedT[t];
      massT[0] += massT[t];
    }
      
    for (int t=1; t<coefs.size(); t++) *coefs[0] += *coefs[t];

    if (use_mpi) {
      MPI_Allreduce(MPI_IN_PLACE, coefs[0]->data(), coefs[0]->size(),
		    MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(&usedT[0], &used,      1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&massT[0], &totalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } else {
      used = usedT[0];
      totalMass = massT[0];
    }
  }


  //! Retrieve the coefficients 
  CoefClasses::CoefStrPtr FieldBasis::getCoefficients()
  {
    if (dof==2) {
      auto cstr = std::make_shared<CoefClasses::PolarVelStruct>();
      cstr->assign(*coefs[0], lmax, nmax);
      return cstr;
    } else {
      auto cstr = std::make_shared<CoefClasses::SphVelStruct>();
      cstr->assign(*coefs[0], lmax, nmax);
      return cstr;
    }
  }
  
  std::vector<double>
  FieldBasis::sph_eval(double r, double costh, double phi)
  {
    constexpr std::complex<double> I(0, 1);

    if (dof==2) {
      std::vector<double> ret(3, 0);
      auto p = (*ortho)(r*sqrt(fabs(1.0 - costh*costh)));
      for (int i=0; i<3; i++) {
	for (int m=0; m<=lmax; m++) {
	  for (int n=0; n<nmax; n++) {
	    ret[i] += std::real((*coefs[0])(i, m, n)*exp(I*(phi*m)))*p[n];
	  }
	}
      }
      return ret;
    } else {
      std::vector<double> ret(4, 0);
      auto p = (*ortho)(r);
      for (int i=0; i<4; i++) {
	for (int l=0, L=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++, L++) {
	    for (int n=0; n<nmax; n++) {
	      ret[i] += std::real((*coefs[0])(i, L, n)*exp(I*(phi*m)))*p[n]*
		Ylm01(l, m)*plgndr(l, m, costh);
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
      coef = std::make_shared<CoefClasses::SphVelStruct>();
    else if (dof==2)
      coef = std::make_shared<CoefClasses::PolarVelStruct>();
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
      coefret = std::make_shared<CoefClasses::SphVelStruct>();
    else if (dof==2)
      coefret = std::make_shared<CoefClasses::PolarVelStruct>();
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
  void FieldBasis::addFromArray(Eigen::VectorXd& m, RowMatrixXd& p, bool roundrobin)
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

    if (p.rows() < 10 and p.cols() > p.rows()) {
      std::cout << "Basis::addFromArray: interpreting your "
		<< p.rows() << "X" << p.cols() << " input array as "
		<< p.cols() << "X" << p.rows() << "." << std::endl;

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
    
    double R   = sqrt(x*x + y*y);
    double r   = sqrt(R*R + z*z);
    double phi = atan2(y, x);
    double cth = z/(r + 1.0e-18);
    
    double vr  = (u*x + v*y)/(R + 1.0e-18);
    double vp  = (u*y - v*x)/(R + 1.0e-18);
    
    return {vr, vp};
  }

  std::vector<double> sphVel(double mass,
			     FieldBasis::PS3& pos, FieldBasis::PS3& vel)
  {
    auto [x, y, z] = pos;
    auto [u, v, w] = vel;
    
    double R   = sqrt(x*x + y*y);
    double r   = sqrt(R*R + z*z);
    double phi = atan2(y, x);
    double cth = z/(r + 1.0e-18);
    
    double vr  = (u*x + v*y + w*z)/(r + 1.0e-18);
    double vt  = (u*z*x + v*z*y - w*R)/(R + 1.0e-18)/(r + 1.0-18);
    double vp  = (u*y - v*x)/(R + 1.0e-18);

    return {vr, vt, vp};
  }

  std::vector<double> crtVel(double mass,
			     FieldBasis::PS3& pos, FieldBasis::PS3& vel)
  {
    return {vel[0], vel[1], vel[2]};
  }


  void VelocityBasis::assignFunc()
  {
    fieldLabels.clear();

    // Field labels (force field labels added below)
    //
    fieldLabels.push_back("density");

    if (coordinates == Coord::Cylindrical) {
      fieldLabels.push_back("v_R");
      fieldLabels.push_back("v_p");
      fieldLabels.push_back("v_z");

      fieldFunc = cylVel;
    } else if (coordinates == Coord::Cartesian) {
      fieldLabels.push_back("v_x");
      fieldLabels.push_back("v_y");
      fieldLabels.push_back("v_z");

      fieldFunc = crtVel;
    } else if (coordinates == Coord::None) {
      fieldLabels.push_back("v_x");
      fieldLabels.push_back("v_y");
      fieldLabels.push_back("v_z");

      fieldFunc = crtVel;
    } else {
      fieldLabels.push_back("v_r");
      fieldLabels.push_back("v_t");
      fieldLabels.push_back("v_p");

      fieldFunc = sphVel;
    }
  }

  VelocityBasis::VelocityBasis(const YAML::Node& conf) : FieldBasis(conf)
  {
    assignFunc();
  }

  VelocityBasis::VelocityBasis(const std::string& confstr) : FieldBasis(confstr)
  {
    assignFunc();
  }


}
// END namespace BasisClasses
