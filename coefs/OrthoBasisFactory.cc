#include <algorithm>

#include <YamlCheck.H>
#include <EXPException.H>
#include <OrthoBasisFactory.H>
#include <DiskModels.H>
#include <exputils.H>
#include <gaussQ.H>

#ifdef HAVE_FE_ENABLE
#include <cfenv>
#endif

extern double Ylm01(int ll, int mm);
extern double plgndr(int l, int m, double x);

namespace OrthoBasisClasses
{
  std::map<OrthoBasis::Coord, std::string> OrthoBasis::coordLabels =
    { {OrthoBasis::Coord::Spherical,   "Spherical"},
      {OrthoBasis::Coord::Cylindrical, "Cylindrical"},
      {OrthoBasis::Coord::Cartesian,   "Cartesian"},
      {OrthoBasis::Coord::None,        "None"} };

  OrthoBasis::OrthoBasis(const YAML::Node& CONF)
  {
    // Copy the YAML config
    //
    node = CONF;

    // Complete the initialization
    //
    initialize();
  }

  OrthoBasis::OrthoBasis(const std::string& confstr)
  {
    try {
      // Read the YAML from a string
      //
      node = YAML::Load(confstr);
    }
    catch (const std::runtime_error& error) {
      std::cout << "Basis constructor: found a problem in the YAML config"
		<< std::endl;
      throw;
    }

    // Complete the initialization
    //
    initialize();
  }

  void OrthoBasis::initialize()
  {
#ifdef HAVE_FE_ENABLE
    // Flag invalid FP results only, such as 0/0 or infinity - infinity
    // or sqrt(-1).
    //
    // feenableexcept(FE_INVALID);
#endif
  
    // Check whether MPI is initialized
    //
    int flag;
    MPI_Initialized(&flag);
    if (flag) use_mpi = true;
    else      use_mpi = false;
    
    // Fall back sanity (works for me but this needs to be fixed
    // generally)
    //
    if (use_mpi) {
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    }

    // Parameters for force
    //
    try {
      // First get the id name; we know it exists because it has been
      // parsed by the factory
      name = node["id"].as<std::string>();
      // Then . . . 
      conf = node["parameters"];
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing force 'parameters' for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << node                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Basis: error parsing YAML");
    }
    
    // Set coefficient center to zero by default
    //
    coefctr = {0.0, 0.0, 0.0};

    // Default null coordinate type
    //
    coordinates = Coord::None;
  }
  
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

  OrthoBasis::Coord OrthoBasis::parseFieldType(std::string coord_type)
  {
    // Find coordinate type
    //
    std::transform(coord_type.begin(), coord_type.end(), coord_type.begin(),
    [](unsigned char c){ return std::tolower(c); });

    Coord ctype;

    if (coord_type.find("cyl") != std::string::npos) {
      ctype = Coord::Cylindrical;
    } else if (coord_type.find("cart") != std::string::npos) {
      ctype = Coord::Cartesian;
    } else if (coord_type.find("none") != std::string::npos) {
      ctype = Coord::None;
    } else {
      ctype = Coord::Spherical;
    }

    return ctype;
  }

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


  // Derived construcgtor for kinematic bases
  //
  VelocityBasis::VelocityBasis(const YAML::Node& CONF) : OrthoBasis(CONF)
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
	return exp(-r/ascl) * 0.5*(1.0 + std::erf((rmax - 5.0*delta - r)/delta)) / ascl;
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
    /*
    // Check if file exists
    //
    if (std::filesystem::exists(outfile)) {
      ExtendH5Coefs();
    }
    // Otherwise, extend the existing HDF5 file
    //
    else {
      WriteH5Coefs();
    }
    */
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


}
// END namespace OrthoBasisClasses
