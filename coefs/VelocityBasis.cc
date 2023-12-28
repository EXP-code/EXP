#include <algorithm>

#include <VelocityBasis.H>
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
  VelocityBasis::valid_keys = {
    "filename",
    "nint",
    "nintsub",
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

  std::vector<std::string>
  VelocityBasis::getFieldLabels(const Coord ctype)
  {
    // Field labels (force field labels added below)
    //
    std::vector<std::string> labels;

    if (ctype == Coord::Cylindrical) {
      labels.push_back("v_R");
      labels.push_back("v_p");
      labels.push_back("v_z");
    } else if (ctype == Coord::Cartesian) {
      labels.push_back("v_x");
      labels.push_back("v_y");
      labels.push_back("v_z");
    } else if (ctype == Coord::None) {
      // No forces
    } else {
      labels.push_back("v_r");
      labels.push_back("v_t");
      labels.push_back("v_p");
    }

    return labels;
  }


  // Deploy instance based on configuration file 
  //
  void VelocityBasis::configure()
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
      sout << "VelocityBasis: found type <" << model << ">.  Must be one of ";
      for (auto v : modelTypes) sout << " " << v;
      throw std::runtime_error(sout.str());
    }
    
    // Check dof value
    //
    if (dof!=2 or dof!=3) {
      std::ostringstream sout;
      sout << "VelocityBasis: found " << dof << " for dof.  Must be 2 or 3.";
      throw std::runtime_error(sout.str());
    }
    
    // Allocate storage for coefficients
    //
    store.resize(omp_get_max_threads());
    coefs.resize(omp_get_max_threads());
    massT.resize(omp_get_max_threads());
    usedT.resize(omp_get_max_threads());
    
    for (int t=0; t<omp_get_max_threads(); t++) {
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
      throw InternalError("VelocityBasis: model logic failure?! "
			  "You should not be here...",
			  __FILE__, __LINE__);
    }
    
    // Generate the orthogonal function instance
    //
    ortho = std::make_shared<OrthoFunction>(nmax, densfunc, rmin, rmax, scale);
  }

  void VelocityBasis::initialize()
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
      if (myid==0) std::cout << "Error parsing parameters in VelocityBasis: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      MPI_Finalize();
      exit(-1);
    }
  }

  void VelocityBasis::load_coefs(CoefClasses::CoefStrPtr cf, double time)
  {
    coefstr = cf;		// Assign the coeffiicient container
    coefstr->time = time;	// and the time

    // Copy the data structure
    if (dof==2) {
      coefs[0] = dynamic_pointer_cast<CoefClasses::SphVelStruct>(cf)->coefs;
    } else {
      coefs[0] = dynamic_pointer_cast<CoefClasses::SphVelStruct>(cf)->coefs;
    }
  }

  void VelocityBasis::reset_coefs()
  {
    used = 0;
    totalMass = 0.0;
    for (auto & v : usedT) v = 0;
    for (auto & v : massT) v = 0.0;
    for (auto & v : store) v.setZero();
  }

  void VelocityBasis::set_coefs(CoefClasses::CoefStrPtr c)
  {
    if (dof==2) {
      auto p = dynamic_cast<CoefClasses::PolarVelStruct*>(c.get());
      coefs[0] = p->coefs;
    } else {
      auto p = dynamic_cast<CoefClasses::SphVelStruct*>(c.get());
      coefs[0] = p->coefs;
    }
  }

  void VelocityBasis::accumulate(double mass,
				 double x, double y, double z,
				 double u, double v, double w)
  {
    constexpr std::complex<double> I(0, 1);
    int tid = omp_get_thread_num();

    double R   = sqrt(x*x + y*y);
    double r   = sqrt(R*R + z*z);
    double phi = atan2(y, x);
    double cth = z/(r + 1.0e-18);
    
    if (dof==2) {
      
      double vr  = (u*x + v*y)/(R + 1.0e-18);
      double vp  = (u*y - v*x)/(R + 1.0e-18);
    
      auto p = (*ortho)(R);
    
      for (int m=0; m<=lmax; m++) {
	std::complex<double> P = std::exp(I*(phi*m));
	for (int n=0; n<nmax; n++) {
	  (*coefs[tid])(0, m, n) += mass*P*p(n);
	  (*coefs[tid])(1, m, n) += mass*P*p(n)*vr;
	  (*coefs[tid])(2, m, n) += mass*P*p(n)*vp;
	}
      }	 
      
    } else {
    
      double vr  = (u*x + v*y + w*z)/(r + 1.0e-18);
      double vt  = (u*z*x + v*z*y - w*R)/(R + 1.0e-18)/(r + 1.0-18);
      double vp  = (u*y - v*x)/(R + 1.0e-18);
    
      auto p = (*ortho)(r);
      
      for (int l=0, k=0; l<=lmax; l++) {
	for (int m=0; m<=lmax; m++, k++) {
	  std::complex<double> P =
	    std::exp(I*(phi*m))*Ylm01(l, m)*plgndr(l, m, cth);
	
	  for (int n=0; n<nmax; n++) {
	    (*coefs[tid])(0, m, n) += mass*P*p(n);
	    (*coefs[tid])(1, m, n) += mass*P*p(n)*vr;
	    (*coefs[tid])(2, m, n) += mass*P*p(n)*vt;
	    (*coefs[tid])(3, m, n) += mass*P*p(n)*vp;
	  }
	}
      }
    }
  }

  void VelocityBasis::make_coefs()
  {
    // Sum over threads
    //
    for (int t=1; t<coefs.size(); t++) {
      *coefs[0] += *coefs[t];
      usedT[0] += usedT[t];
      massT[0] += massT[t];
    }
      
    for (int t=1; t<coefs.size(); t++) *coefs[0] += *coefs[t];

    MPI_Allreduce(MPI_IN_PLACE, coefs[0]->data(), coefs[0]->size(),
		  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(&usedT[0], &used,      1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&massT[0], &totalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }


  //! Retrieve the coefficients 
  CoefClasses::CoefStrPtr VelocityBasis::getCoefficients()
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
  VelocityBasis::sph_eval(double r, double costh, double phi)
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
      for (int i=0; i<3; i++) {
	for (int l=0, L=0; l<=lmax; l++) {
	  for (int m=0; m<=lmax; m++, L++) {
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
  VelocityBasis::cyl_eval(double R, double z, double phi)
  {
    double r = sqrt(R*R + z*z);
    double costh = z/(r + 1.0e-18);

    return sph_eval(r, costh, phi);
  }

  std::vector<double>
  VelocityBasis::crt_eval(double x, double y, double z)
  {
    double r = sqrt(x*x + y*y + z*z);
    double costh = z/(r + 1.0e-18), phi = atan2(y, x);
    
    return sph_eval(r, costh, phi);
  }

  VelocityBasis::BasisArray
  VelocityBasis::getBasis(double logxmin, double logxmax, int numgrid)
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
}
// END namespace OrthoBasisClasses
