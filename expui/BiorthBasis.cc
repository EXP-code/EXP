#include <algorithm>
#include <set>

#include <YamlCheck.H>
#include <EXPException.H>
#include <BiorthBasis.H>
#include <DiskModels.H>
#include <config_exp.h>
#include <exputils.H>
#include <gaussQ.H>

#ifdef HAVE_FE_ENABLE
#include <cfenv>
#endif

namespace BasisClasses
{
  const std::set<std::string>
  Spherical::valid_keys = {
    "rmapping",
    "cmap",
    "Lmax",
    "dof",
    "npca",
    "npca0",
    "pcavar",
    "pcadiag",
    "pcavtk",
    "pcaeof",
    "subsamp",
    "sampT",
    "snr",
    "samplesz",
    "vtkfreq",
    "tksmooth",
    "tkcum",
    "tk_type",
    "nmax",
    "scale",
    "rmin",
    "rmax",
    "numr",
    "nums",
    "N1",
    "N2",
    "NO_L0",
    "NO_L1",
    "EVEN_L",
    "EVEN_M",
    "M0_ONLY",
    "NOISE",
    "noiseN",
    "noise_model_file",
    "seedN",
    "ssfrac",
    "playback",
    "coefCompute",
    "coefMaster",
    "diverge",
    "dfac",
    "dtime",
    "logr",
    "plummer",
    "self_consistent",
    "cachename",
    "modelname",
    "pyname",
    "rnum"
  };

  std::vector<std::string> BiorthBasis::getFieldLabels(const Coord ctype)
  {
    // Field labels (force field labels added below)
    //
    std::vector<std::string> labels =
      {"dens m=0", "dens m>0", "dens",
       "potl m=0", "potl m>0", "potl"};

    if (ctype == Coord::Cylindrical) {
      labels.push_back("rad force");
      labels.push_back("ver force");
      labels.push_back("azi force");
      if (midplane) labels.push_back("midplane");
    } else if (ctype == Coord::Cartesian) {
      labels.push_back("x force");
      labels.push_back("y force");
      labels.push_back("z force");
    } else if (ctype == Coord::None) {
      // No forces
    } else {
      labels.push_back("rad force");
      labels.push_back("mer force");
      labels.push_back("azi force");
    }

    return labels;
  }

  std::shared_ptr<BiorthBasis> BiorthBasis::factory_string(const std::string& conf)
  {
    YAML::Node node;
    
    try {
      // Read the YAML from a string
      //
      node = YAML::Load(conf);
    }
    catch (const std::runtime_error& error) {
      std::cout << "BiorthBasis::factory constructor: found a problem in the YAML config"
		<< std::endl;
      throw;
    }

    // Complete the initialization
    //
    return factory_initialize(node);
  }

  std::shared_ptr<BiorthBasis> BiorthBasis::factory(const YAML::Node& conf)
  {
    return factory_initialize(conf);
  }

  std::shared_ptr<BiorthBasis> BiorthBasis::factory_initialize(const YAML::Node& conf)
  {
    std::shared_ptr<BiorthBasis> basis;
    std::string name;
    
    // Load parameters from YAML configuration node
    try {
      name = conf["id"].as<std::string>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing force id in BiorthBasis::factory"
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("BiorthBasis::factory: error parsing YAML");
    }
    
    try {
      if ( !name.compare("sphereSL") ) {
	basis = std::make_shared<SphericalSL>(conf);
      }
      else if ( !name.compare("bessel") ) {
	basis = std::make_shared<Bessel>(conf);
      }
      else if ( !name.compare("cylinder") ) {
	basis = std::make_shared<Cylindrical>(conf);
      }
      else if ( !name.compare("flatdisk") ) {
	basis = std::make_shared<FlatDisk>(conf);
      }
      else if ( !name.compare("CBDisk") ) {
	basis = std::make_shared<CBDisk>(conf);
      }
      else if ( !name.compare("slabSL") ) {
	basis = std::make_shared<Slab>(conf);
      }
      else if ( !name.compare("cube") ) {
	basis = std::make_shared<Cube>(conf);
      }
      else {
	std::string msg("I don't know about the basis named: ");
	msg += name;
	msg += ". Known types are currently 'sphereSL', 'bessel', 'cylinder', 'flatdisk', 'CBDisk', 'slabSL', and 'cube'";
	throw std::runtime_error(msg);
      }
    }
    catch (std::exception& e) {
      std::cout << "Error in BiorthBasis::factory constructor: " << e.what() << std::endl;
      throw;			// Rethrow the exception?
    }
    
    return basis;
  }


  Spherical::Spherical(const YAML::Node& CONF, const std::string& forceID) :
    BiorthBasis(CONF, forceID)
  {
    initialize();
  }

  Spherical::Spherical(const std::string& confstr, const std::string& forceID) :
    BiorthBasis(confstr, forceID)
  {
    initialize();
  }

  SphericalSL::SphericalSL(const YAML::Node& CONF) : Spherical(CONF, "sphereSL")
  {
    initialize();
  }

  SphericalSL::SphericalSL(const std::string& confstr) : Spherical(confstr, "sphereSL")
  {
    initialize();
  }

  Bessel::Bessel(const YAML::Node& CONF) :  Spherical(CONF, "Bessel")
  {
    initialize();
  }

  Bessel::Bessel(const std::string& confstr) : Spherical(confstr, "Bessel")
  {
    initialize();
  }

  void Spherical::initialize()
  {

    // Assign some defaults
    //
    cmap       = 1;
    lmax       = 6;
    nmax       = 18;
    
    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::Spherical", "parameter", unmatched, __FILE__, __LINE__);
    
    try {
      if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
      if (conf["Lmax"])      lmax       = conf["Lmax"].as<int>();
      if (conf["nmax"])      nmax       = conf["nmax"].as<int>();
      
      if (conf["rmapping"]) 
	rmap = conf["rmapping"].as<double>();
      else
	rmap = 1.0;
      
      if (conf["scale"]) 
	scale = conf["scale"].as<double>();
      else
	scale = 1.0;
      
      if (conf["rmin"]) 
	rmin = conf["rmin"].as<double>();
      else
	rmin = 0.0;
      
      if (conf["rmax"]) 
	rmax = conf["rmax"].as<double>();
      else
	rmax = std::numeric_limits<double>::max();
      
      if (conf["numr"])
	numr = conf["numr"].as<int>();
      else
	numr = 800;
      
      N1 = 0;
      N2 = std::numeric_limits<int>::max();
      NO_L0 = NO_L1 = EVEN_L = EVEN_M = M0_only = false;
      
      if (conf["N1"]   )     N1        = conf["N1"].as<bool>();
      if (conf["N2"]   )     N2        = conf["N2"].as<bool>();
      if (conf["NO_L0"])     NO_L0     = conf["NO_L0"].as<bool>();
      if (conf["NO_L1"])     NO_L1     = conf["NO_L1"].as<bool>();
      if (conf["EVEN_L"])    EVEN_L    = conf["EVEN_L"].as<bool>();
      if (conf["EVEN_M"])    EVEN_M    = conf["EVEN_M"].as<bool>();
      if (conf["M0_ONLY"])   M0_only   = conf["M0_ONLY"].as<bool>();
      if (conf["pcavar"])    pcavar    = conf["pcavar"].as<bool>();
      if (conf["sampT"])     sampT     = conf["sampT"].as<int>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Spherical: error parsing YAML");
    }

    // Number of possible threads
    int nthrds = omp_get_max_threads();
    
    potd.resize(nthrds);
    dpot.resize(nthrds);
    dpt2.resize(nthrds);
    dend.resize(nthrds);
    
    legs  .resize(nthrds);
    dlegs .resize(nthrds);
    d2legs.resize(nthrds);

    for (auto & v : potd) v.resize(lmax+1, nmax);
    for (auto & v : dpot) v.resize(lmax+1, nmax);
    for (auto & v : dpt2) v.resize(lmax+1, nmax);
    for (auto & v : dend) v.resize(lmax+1, nmax);

    // Wasteful but simple factorial table.  Could be done with a
    // triangular indexing scheme...
    //
    for (auto & v : legs  ) v.resize(lmax+1, lmax+1);
    for (auto & v : dlegs ) v.resize(lmax+1, lmax+1);
    for (auto & v : d2legs) v.resize(lmax+1, lmax+1);

    // Allocate coefficient storage; stores real and imaginary parts as reals
    //
    expcoef.resize((lmax+1)*(lmax+1), nmax);
    expcoef.setZero();
      
    work.resize(nmax);
      
    // Wasteful but simple factorial table.  Could be done with a
    // triangular indexing scheme...
    //
    factorial.resize(lmax+1, lmax+1);
      
    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	factorial(l, m) = sqrt( (0.5*l+0.25)/M_PI * 
				exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );
	if (m != 0) factorial(l, m) *= M_SQRT2;
      }
    }

    used = 0;

    // Initialize covariance
    //
    if (pcavar) init_covariance();

    // Set spherical coordindates
    //
    coordinates = Coord::Spherical;
  }
  
  void Spherical::init_covariance()
  {
    if (pcavar) {
      // Triangular for l,m>=0
      int Ltot = (lmax+1)*(lmax+2)/2;
	
      meanV.resize(sampT);
      for (auto& v : meanV) {
	v.resize(Ltot);
	for (auto& vec : v) vec.resize(nmax);
      }

      covrV.resize(sampT);
      for (auto& v : covrV) {
	v.resize(Ltot);
	for (auto& mat : v) mat.resize(nmax, nmax);
      }

      sampleCounts.resize(sampT);
      sampleMasses.resize(sampT);
      
      zero_covariance();
    }
  }

  void Spherical::zero_covariance()
  {
    for (int T=0; T<sampT; T++) {
      for (auto& v : meanV[T]) v.setZero();
      for (auto& v : covrV[T]) v.setZero();
    }

    sampleCounts.setZero();
    sampleMasses.setZero();
  }

  void SphericalSL::initialize()
  {
    // Identifier
    //
    BasisID = "sphereSL";

    // Assign some defaults
    //
    model_file = "SLGridSph.model";
    
    // Default cachename, empty by default
    //
    std::string cachename;

    try {
      if (conf["modelname"]) model_file = conf["modelname"].as<std::string>();
      if (conf["cachename"]) cachename  = conf["cachename"].as<std::string>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("SphericalSL: error parsing YAML");
    }
    
    // Check for non-null cache file name.  This must be specified
    // to prevent recomputation and unexpected behavior.
    //
    if (cachename.size() == 0) {
      throw std::runtime_error
	("SphericalSL requires a specified cachename in your YAML config\n"
	 "for consistency with previous invocations and existing coefficient\n"
	 "sets.  Please add explicitly add 'cachename: name' to your config\n"
	 "with new 'name' for creating a basis or an existing 'name' for\n"
	 "reading a previously generated basis cache\n");
    }

    // Set MPI flag in SLGridSph from MPI_Initialized
    SLGridSph::mpi = use_mpi ? 1 : 0;
    
    // Instantiate to get min/max radius from the model
    mod = std::make_shared<SphericalModelTable>(model_file);
    
    // Set rmin to a sane value if not specified
    if (not conf["rmin"] or rmin < mod->get_min_radius()) 
      rmin = mod->get_min_radius();

    // Set rmax to a sane value if not specified
    if (not conf["rmax"] or rmax > mod->get_max_radius()) 
      rmax = mod->get_max_radius()*0.99;
    
    // Finally, make the Sturm-Lioville basis...
    sl = std::make_shared<SLGridSph>
      (model_file, lmax, nmax, numr, rmin, rmax, true, cmap, rmap,
       0, 1, cachename);
    
    // Test basis for consistency
    if (myid==0) orthoTest(200);
  }
  
  void Bessel::initialize()
  {
    // Identifier
    //
    BasisID = "bessel";

    try {
      if (conf["rnum"])
	rnum = conf["rnum"].as<int>();
      else
	rnum = 2000;
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("SphericalSL: error parsing YAML");
    }
    
    // Finally, make the Sturm-Lioville basis...
    bess = std::make_shared<BiorthBess>(rmax, lmax, nmax, rnum);

    // Test basis for consistency
    if (myid==0) orthoTest(200);
  }
  
  void Spherical::reset_coefs(void)
  {
    if (expcoef.rows()>0 && expcoef.cols()>0) expcoef.setZero();
    totalMass = 0.0;
    used = 0;
    if (pcavar) zero_covariance();
  }
  
  
  void Spherical::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    CoefClasses::SphStruct* cf = dynamic_cast<CoefClasses::SphStruct*>(coef.get());

    cf->lmax   = lmax;
    cf->nmax   = nmax;
    cf->scale  = scale;
    cf->time   = time;
    cf->normed = true;

    // Angular storage dimension; triangular number for complex coefs
    int ldim = (lmax+1)*(lmax+2)/2;

    // Allocate the coefficient storage
    cf->store.resize(ldim*nmax);

    // Make the coefficient map
    cf->coefs = std::make_shared<CoefClasses::SphStruct::coefType>
      (cf->store.data(), ldim, nmax);

    for (int l=0, L0=0, L1=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<nmax; n++) {
	  if (m==0)
	    (*cf->coefs)(L0, n) = {expcoef(L1, n), 0.0};
	  else
	    (*cf->coefs)(L0, n) = {expcoef(L1, n), expcoef(L1+1, n)};
	}
	L0 += 1;
	if (m==0) L1 += 1;
	else      L1 += 2;
      }
    }
  }

  void Spherical::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (not dynamic_cast<CoefClasses::SphStruct*>(coef.get()))
      throw std::runtime_error("Spherical::set_coefs: you must pass a CoefClasses::SphStruct");

    // Sanity check on dimensionality
    //
    {
      auto & cf = *dynamic_cast<CoefClasses::SphStruct*>(coef.get());

      int rows = (*cf.coefs).rows();
      int cols = (*cf.coefs).cols();
      int rexp = (lmax+1)*(lmax+2)/2;
      if (rows != rexp or cols != nmax) {
	std::ostringstream sout;
	sout << "Spherical::set_coefs: the basis has (lmax, nmax)=("
	     << lmax << ", " << nmax
	     << ") and the dimensions must be (rows, cols)=("
	     << rexp << ", " << nmax
	     << "). The coef structure has (rows, cols)=("
	     << rows << ", " << cols << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    CoefClasses::SphStruct* cf = dynamic_cast<CoefClasses::SphStruct*>(coef.get());

    // Cache the current coefficient structure
    //
    coefret = coef;

    // Assign internal coefficient table (doubles) from the complex struct
    //
    for (int l=0, L0=0, L1=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<nmax; n++) {
	  if (m==0)
	    expcoef(L1,   n) = (*cf->coefs)(L0, n).real();
	  else {
	    expcoef(L1,   n) = (*cf->coefs)(L0, n).real();
	    expcoef(L1+1, n) = (*cf->coefs)(L0, n).imag();
	  }
	}
	L0 += 1;
	if (m==0) L1 += 1;
	else      L1 += 2;
      }
    }

    // Assign center if need be
    //
    if (cf->ctr.size())
      coefctr = cf->ctr;
    else
      coefctr = {0.0, 0.0, 0.0};
  }

  void Spherical::accumulate(double x, double y, double z, double mass)
  {
    double fac, fac1, fac2, fac4;
    double norm = -4.0*M_PI;
    const double dsmall = 1.0e-20;
    
    // Get thread id
    int tid = omp_get_thread_num();

    //======================
    // Compute coefficients 
    //======================
    
    double r2 = (x*x + y*y + z*z);
    double r = sqrt(r2) + dsmall;
    double costh = z/r;
    double phi = atan2(y,x);
    double rs = r/scale;
    
    if (r < rmin or r > rmax) return;
    
    used++;
    totalMass += mass;
    
    get_pot(potd[tid], rs);
    
    legendre_R(lmax, costh, legs[tid]);
    
    // Sample index for pcavar
    int T = 0;
    if (pcavar) {
      T = used % sampT;
      sampleCounts(T) += 1;
      sampleMasses(T) += mass;
    }
    
    // L loop
    for (int l=0, loffset=0, L=0; l<=lmax; loffset+=(2*l+1), l++) {
      
      Eigen::VectorXd workE;
      int esize = (l+1)*nmax;
      
      // M loop
      for (int m=0, moffset=0, moffE=0; m<=l; m++, L++) {
	
	if (m==0) {
	  fac = factorial(l, m) * legs[tid](l, m);
	  for (int n=0; n<nmax; n++) {
	    fac4 = potd[tid](l, n)*fac;
	    expcoef(loffset+moffset, n) += fac4 * norm * mass;

	    if (pcavar) {
	      meanV[T][L](n) += fac4 * norm * mass;
	      for (int o=0; o<nmax; o++)
		covrV[T][L](n, o)  += fac4 * norm *
		  fac * potd[tid](l, o) * norm * mass;
	    }
	  }
	  
	  moffset++;
	}
	else {
	  fac  = factorial(l, m) * legs[tid](l, m);
	  fac1 = fac*cos(phi*m);
	  fac2 = fac*sin(phi*m);
	  for (int n=0; n<nmax; n++) {
	    fac4 = potd[tid](l, n);
	    expcoef(loffset+moffset  , n) += fac1 * fac4 * norm * mass;
	    expcoef(loffset+moffset+1, n) += fac2 * fac4 * norm * mass;

	    if (pcavar) {
	      meanV[T][L](n) += fac * fac4 * norm * mass;
	      for (int o=0; o<nmax; o++)
		covrV[T][L](n, o)  += fac * fac4 * norm *
		  fac * potd[tid](l, o) * norm * mass;
	    }
	  }
	  
	  moffset+=2;
	}
      }
    }
    
  }
  
  void Spherical::make_coefs()
  {
    // MPI reduction of coefficients
    if (use_mpi) {
      
      // Square of total number of angular coefficients in real form
      int Ltot = (lmax+1)*(lmax+1);

      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      for (int l=0; l<Ltot; l++) {
	work = expcoef.row(l);
	MPI_Allreduce(MPI_IN_PLACE, work.data(), nmax, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	expcoef.row(l) = work;
      }
      
      if (pcavar) {

	// Triangular number of angular coefficients in complex form
	int ltot = (lmax+1)*(lmax+2)/2;

	MPI_Allreduce(MPI_IN_PLACE, sampleCounts.data(), sampleCounts.size(),
		      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	MPI_Allreduce(MPI_IN_PLACE, sampleMasses.data(), sampleMasses.size(),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for (int T=0; T<sampT; T++) {

	  for (int l=0; l<ltot; l++) {
	    
	    MPI_Allreduce(MPI_IN_PLACE, meanV[T][l].data(), meanV[T][l].size(),
			  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    MPI_Allreduce(MPI_IN_PLACE, covrV[T][l].data(), covrV[T][l].size(),
			  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  }
	}
      }
    }
  }
  
  std::vector<double>
  Spherical::sph_eval(double r, double costh, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    double fac1 = factorial(0, 0);
    
    get_dens (dend[tid], r/scale);
    get_pot  (potd[tid], r/scale);
    get_force(dpot[tid], r/scale);
    
    legendre_R(lmax, costh, legs[tid], dlegs[tid]);
    
    double den0, pot0, potr;

    if (NO_L0) {
      den0 = 0.0;
      pot0 = 0.0;
      potr = 0.0;
    } else {
      den0 = fac1 * expcoef.row(0).dot(dend[tid].row(0));
      pot0 = fac1 * expcoef.row(0).dot(potd[tid].row(0));
      potr = fac1 * expcoef.row(0).dot(dpot[tid].row(0));
    }

    double den1 = 0.0;
    double pot1 = 0.0;
    double pott = 0.0;
    double potp = 0.0;
    
    // L loop
    for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
      
      // Check for even l
      if (EVEN_L and l%2) continue;
      
      // No l=1
      if (NO_L1 and l==1) continue;
      
      // M loop
      for (int m=0, moffset=0; m<=l; m++) {
	
	if (M0_only and m) continue;
	if (EVEN_M and m%2) continue;
	
	fac1 = factorial(l, m);
	if (m==0) {
	  double sumR=0.0, sumP=0.0, sumD=0.0;
	  for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	    sumR += expcoef(loffset+moffset, n) * dend[tid](l, n);
	    sumP += expcoef(loffset+moffset, n) * potd[tid](l, n);
	    sumD += expcoef(loffset+moffset, n) * dpot[tid](l, n);
	  }
	  
	  den1 += fac1*legs[tid] (l, m) * sumR;
	  pot1 += fac1*legs[tid] (l, m) * sumP;
	  potr += fac1*legs[tid] (l, m) * sumD;
	  pott += fac1*dlegs[tid](l, m) * sumP;
	  
	  moffset++;
	}
	else {
	  double cosm = cos(phi*m);
	  double sinm = sin(phi*m);
	  
	  double sumR0=0.0, sumP0=0.0, sumD0=0.0;
	  double sumR1=0.0, sumP1=0.0, sumD1=0.0;
	  for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	    sumR0 += expcoef(loffset+moffset+0, n) * dend[tid](l, n);
	    sumP0 += expcoef(loffset+moffset+0, n) * potd[tid](l, n);
	    sumD0 += expcoef(loffset+moffset+0, n) * dpot[tid](l, n);
	    sumR1 += expcoef(loffset+moffset+1, n) * dend[tid](l, n);
	    sumP1 += expcoef(loffset+moffset+1, n) * potd[tid](l, n);
	    sumD1 += expcoef(loffset+moffset+1, n) * dpot[tid](l, n);
	  }
	  
	  den1 += fac1 * legs[tid] (l, m) *  ( sumR0*cosm + sumR1*sinm );
	  pot1 += fac1 * legs[tid] (l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potr += fac1 * legs[tid] (l, m) *  ( sumD0*cosm + sumD1*sinm );
	  pott += fac1 * dlegs[tid](l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potp += fac1 * legs[tid] (l, m) *  (-sumP0*sinm + sumP1*cosm ) * m;
	  
	  moffset +=2;
	}
      }
    }
    
    double densfac = 1.0/(scale*scale*scale) * 0.25/M_PI;
    double potlfac = 1.0/scale;

    return
      {den0 * densfac,		// 0
       den1 * densfac,		// 1
       (den0 + den1) * densfac,	// 2
       pot0 * potlfac,		// 3
       pot1 * potlfac,		// 4
       (pot0 + pot1) * potlfac,	// 5
       potr * (-potlfac)/scale,	// 6
       pott * (-potlfac),	// 7
       potp * (-potlfac)};	// 8
    //         ^
    //         |
    // Return force not potential gradient
  }

  void Spherical::computeAccel(double x, double y, double z,
			       Eigen::Ref<Eigen::Vector3d> acc)
  {
    // Get polar coordinates
    double R     = sqrt(x*x + y*y);
    double r     = sqrt(R*R + z*z);
    double costh = z/r;
    double sinth = R/r;
    double phi   = atan2(y, x);

    // Get thread id
    int tid = omp_get_thread_num();

    
    double fac1 = factorial(0, 0);
    
    get_pot  (potd[tid], r/scale);
    get_force(dpot[tid], r/scale);
    
    legendre_R(lmax, costh, legs[tid], dlegs[tid]);
    
    double pot0, potr;

    if (NO_L0) {
      pot0 = 0.0;
      potr = 0.0;
    } else {
      pot0 = fac1 * expcoef.row(0).dot(potd[tid].row(0));
      potr = fac1 * expcoef.row(0).dot(dpot[tid].row(0));
    }

    double den1 = 0.0;
    double pot1 = 0.0;
    double pott = 0.0;
    double potp = 0.0;
    
    // L loop
    for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
      
      // Check for even l
      if (EVEN_L and l%2) continue;
      
      // No l=1
      if (NO_L1 and l==1) continue;
      
      // M loop
      for (int m=0, moffset=0; m<=l; m++) {
	
	if (M0_only and m) continue;
	if (EVEN_M and m%2) continue;
	
	fac1 = factorial(l, m);
	if (m==0) {
	  double sumR=0.0, sumP=0.0, sumD=0.0;
	  for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	    sumR += expcoef(loffset+moffset, n) * dend[tid](l, n);
	    sumP += expcoef(loffset+moffset, n) * potd[tid](l, n);
	    sumD += expcoef(loffset+moffset, n) * dpot[tid](l, n);
	  }
	  
	  pot1 += fac1*legs[tid] (l, m) * sumP;
	  potr += fac1*legs[tid] (l, m) * sumD;
	  pott += fac1*dlegs[tid](l, m) * sumP;
	  
	  moffset++;
	}
	else {
	  double cosm = cos(phi*m);
	  double sinm = sin(phi*m);
	  
	  double sumR0=0.0, sumP0=0.0, sumD0=0.0;
	  double sumR1=0.0, sumP1=0.0, sumD1=0.0;
	  for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	    sumP0 += expcoef(loffset+moffset+0, n) * potd[tid](l, n);
	    sumD0 += expcoef(loffset+moffset+0, n) * dpot[tid](l, n);
	    sumP1 += expcoef(loffset+moffset+1, n) * potd[tid](l, n);
	    sumD1 += expcoef(loffset+moffset+1, n) * dpot[tid](l, n);
	  }
	  
	  pot1 += fac1 * legs[tid] (l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potr += fac1 * legs[tid] (l, m) *  ( sumD0*cosm + sumD1*sinm );
	  pott += fac1 * dlegs[tid](l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potp += fac1 * legs[tid] (l, m) *  (-sumP0*sinm + sumP1*cosm ) * m;
	  
	  moffset +=2;
	}
      }
    }
    
    double potlfac = 1.0/scale;

    potr *= (-potlfac)/scale;
    pott *= (-potlfac);
    potp *= (-potlfac);

    double potR = potr*sinth + pott*costh;
    double potz = potr*costh - pott*sinth;
  
    double tpotx = potR*x/R - potp*y/R ;
    double tpoty = potR*y/R + potp*x/R ;

    // Return force not potential gradient
    //
    acc << tpotx, tpoty, potz;
  }


  std::vector<double>
  Spherical::cyl_eval(double R, double z, double phi)
  {
    double r = sqrt(R*R + z*z) + 1.0e-18;
    double costh = z/r, sinth = R/r;
    
    auto v = sph_eval(r, costh, phi);
    
    double potR = v[6]*sinth + v[7]*costh;
    double potz = v[6]*costh - v[7]*sinth;

    return {v[0], v[1], v[2], v[3], v[4], v[5], potR, potz, v[8]};
  }


  // Evaluate in cartesian coordinates
  std::vector<double>
  Spherical::crt_eval
  (double x, double y, double z)
  {
    double R = sqrt(x*x + y*y) + 1.0e-18;
    double phi = atan2(y, x);

    auto v = cyl_eval(R, z, phi);

    double tpotx = v[6]*x/R - v[8]*y/R ;
    double tpoty = v[6]*y/R + v[8]*x/R ;
    
    return {v[0], v[1], v[2], v[3], v[4], v[5], tpotx, tpoty, v[7]};
  }
  
  Spherical::BasisArray SphericalSL::getBasis
  (double logxmin, double logxmax, int numgrid)
  {
    // Assing return storage
    BasisArray ret (lmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numgrid); // Potential
	u["density"  ].resize(numgrid); // Density
	u["rforce"   ].resize(numgrid); // Radial force
      }
    }

    // Radial grid spacing
    double dx = (logxmax - logxmin)/(numgrid-1);

    // Basis storage
    Eigen::MatrixXd tabpot, tabden, tabfrc;

    for (int i=0; i<numgrid; i++) {
      get_pot  (tabpot, pow(10.0, logxmin + dx*i));
      get_dens (tabden, pow(10.0, logxmin + dx*i));
      get_force(tabfrc, pow(10.0, logxmin + dx*i));
      for (int l=0; l<=lmax; l++) {
	for (int n=0; n<nmax;n++){
	  ret[l][n]["potential"](i) = tabpot(l, n);
	  ret[l][n]["density"  ](i) = tabden(l, n);
	  ret[l][n]["rforce"   ](i) = tabfrc(l, n);
	}
      }
    }
    
    return ret;
  }

  Spherical::BasisArray Spherical::getBasis
  (double rmin, double rmax, int numgrid)
  {
    // Assing return storage
    BasisArray ret (lmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numgrid); // Potential
	u["density"  ].resize(numgrid); // Density
	u["rforce"   ].resize(numgrid); // Radial force
      }
    }

    // Radial grid spacing
    double dr = (rmax - rmin)/(numgrid-1);

    // Basis storage
    Eigen::MatrixXd tabpot, tabden, tabfrc;

    for (int i=0; i<numgrid; i++) {
      get_pot  (tabpot, rmin + dr*i);
      get_dens (tabden, rmin + dr*i);
      get_force(tabfrc, rmin + dr*i);
      for (int l=0; l<=lmax; l++) {
	for (int n=0; n<nmax;n++){
	  ret[l][n]["potential"](i) = tabpot(l, n);
	  ret[l][n]["density"  ](i) = tabden(l, n);
	  ret[l][n]["rforce"   ](i) = tabfrc(l, n);
	}
      }
    }
    
    return ret;
  }

  
  /** Return a vector of tuples of basis functions and the covariance
      matrix for subsamples of particles */
  std::vector<std::vector<BiorthBasis::CoefCovarType>>
  Spherical::getCoefCovariance()
  {
    std::vector<std::vector<BiorthBasis::CoefCovarType>> ret;
   if (pcavar) {
     ret.resize(sampT);
     for (int T=0; T<sampT; T++) {
       int ltot = (lmax+1)*(lmax+2)/2;
       ret[T].resize(ltot);
       for (int l=0; l<ltot; l++) {
	 std::get<0>(ret[T][l]) = meanV[T][l];
	 std::get<1>(ret[T][l]) = covrV[T][l];
       }
     }
   }
    
    return ret;
  }


#define MINEPS 1.0e-10
  
  void BiorthBasis::legendre_R(int lmax, double x, Eigen::MatrixXd& p)
  {
    double fact, somx2, pll, pl1, pl2;
    
    p(0, 0) = pll = 1.0;
    if (lmax > 0) {
      somx2 = sqrt( (1.0 - x)*(1.0 + x) );
      fact = 1.0;
      for (int m=1; m<=lmax; m++) {
	pll *= -fact*somx2;
	p(m, m) = pll;
	if (std::isnan(p(m, m))) {
	  std::ostringstream sout;
	  sout << "BiorthBasis::legendre_R NaN: p[" << m << "][" << m << "]: pll=" << pll;
	  throw std::runtime_error(sout.str());
	}
	fact += 2.0;
      }
    }
    
    for (int m=0; m<lmax; m++) {
      pl2 = p(m, m);
      p(m+1, m) = pl1 = x*(2*m+1)*pl2;
      for (int l=m+2; l<=lmax; l++) {
	p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
	if (std::isnan(p(l, m))) {
	  std::ostringstream sout;
	  sout << "BiorthBasis::legendre_R NaN: p[" << l << "][" << m << "]: pll=" << pll;
	  throw std::runtime_error(sout.str());
	}

	pl2 = pl1;
	pl1 = pll;
      }
    }
    
      if (std::isnan(x)) {
	std::ostringstream sout;
	sout << "BiorthBasis::legendre_R NaN: x";
	throw std::runtime_error(sout.str());
      }
    for (int l=0; l<=lmax; l++)
      for (int m=0; m<=l; m++)
	if (std::isnan(p(l, m))) {
	  std::ostringstream sout;
	  sout << "BiorthBasis::legendre_R NaN: p[" << l << "][" << m << "] lmax=" << lmax;
	  throw std::runtime_error(sout.str());
	}
  }
  
  void BiorthBasis::legendre_R(int lmax, double x, Eigen::MatrixXd& p,
			      Eigen::MatrixXd &dp)
  {
    double fact, somx2, pll, pl1, pl2;
    
    p(0, 0) = pll = 1.0;
    if (lmax > 0) {
      somx2 = sqrt( (1.0 - x)*(1.0 + x) );
      fact = 1.0;
      for (int m=1; m<=lmax; m++) {
	pll *= -fact*somx2;
	p(m, m) = pll;
	fact += 2.0;
      }
    }
    
    for (int m=0; m<lmax; m++) {
      pl2 = p(m, m);
      p(m+1, m) = pl1 = x*(2*m+1)*pl2;
      for (int l=m+2; l<=lmax; l++) {
	p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
	pl2 = pl1;
	pl1 = pll;
      }
    }
    
    if (1.0-fabs(x) < MINEPS) {
      if (x>0) x =   1.0 - MINEPS;
      else     x = -(1.0 - MINEPS);
    }
    
    somx2 = 1.0/(x*x - 1.0);
    dp(0, 0) = 0.0;
    for (int l=1; l<=lmax; l++) {
      for (int m=0; m<l; m++)
	dp(l, m) = somx2*(x*l*p(l, m) - (l+m)*p(l-1, m));
      dp(l, l) = somx2*x*l*p(l, l);
    }
  }
  
  void BiorthBasis::legendre_R(int lmax, double x, Eigen::MatrixXd &p,
			      Eigen::MatrixXd &dp, Eigen::MatrixXd& d2p)
  {
    double fact, somx2, pll, pl1, pl2;
    
    p(0, 0) = pll = 1.0;
    if (lmax > 0) {
      somx2 = sqrt( (1.0 - x)*(1.0 + x) );
      fact = 1.0;
      for (int m=1; m<=lmax; m++) {
	pll *= -fact*somx2;
	p(m, m) = pll;
	fact += 2.0;
      }
    }
    
    for (int m=0; m<lmax; m++) {
      pl2 = p(m, m);
      p(m+1, m) = pl1 = x*(2*m+1)*pl2;
      for (int l=m+2; l<=lmax; l++) {
	p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
	pl2 = pl1;
	pl1 = pll;
      }
    }
    
    if (1.0-fabs(x) < MINEPS) {
      if (x>0) x =   1.0 - MINEPS;
      else     x = -(1.0 - MINEPS);
    }
    
    somx2 = 1.0/(x*x - 1.0);
    dp(0, 0) = 0.0;
    if (lmax) {
      for (int l=1; l<=lmax; l++) {
	for (int m=0; m<l; m++)
	  dp(l, m) = somx2*(x*l*p(l, m) - (l+m)*p(l-1, m));
	dp(l, l) = somx2*x*l*p(l, l);
      }
    }
    
    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++)
	d2p(l, m) = -somx2*(2.0*x*dp(l, m) - p(l, m)*(somx2*m*m + l*(l+1)));
    }
    
  }
  
  // Disk target types for cylindrical basis construction
  const std::map<std::string, Cylindrical::DiskType> Cylindrical::dtlookup =
    { {"constant",    DiskType::constant},
      {"gaussian",    DiskType::gaussian},
      {"mn",          DiskType::mn},
      {"exponential", DiskType::exponential},
      {"doubleexpon", DiskType::doubleexpon},
      {"diskbulge",   DiskType::diskbulge},
      {"python",      DiskType::python}
    };

  const std::set<std::string>
  Cylindrical::valid_keys = {
    "tk_type",
    "rcylmin",
    "rcylmax",
    "acyl",
    "hcyl",
    "sech2",
    "snr",
    "evcut",
    "nmaxfid",
    "lmaxfid",
    "mmax",
    "mlim",
    "nmax",
    "ncylnx",
    "ncylny",
    "ncylr",
    "ncylodd",
    "ncylrecomp",
    "npca",
    "npca0",
    "nvtk",
    "cachename",
    "oldcache",
    "eof_file",
    "override",
    "samplesz",
    "rnum",
    "pnum",
    "tnum",
    "ashift",
    "expcond",
    "ignore",
    "deproject",
    "logr",
    "pcavar",
    "pcaeof",
    "pcavtk",
    "pcadiag",
    "subsamp",
    "try_cache",
    "density",
    "EVEN_M",
    "cmap",
    "cmapr",
    "cmapz",
    "aratio",
    "hratio",
    "dweight",
    "Mfac",
    "HERNA",
    "rwidth",
    "rfactor",
    "rtrunc",
    "rpow",
    "mtype",
    "dtype",
    "vflag",
    "self_consistent",
    "playback",
    "coefCompute",
    "coefMaster",
    "pyname"
  };

  Cylindrical::Cylindrical(const YAML::Node& CONF) :
    BiorthBasis(CONF, "cylinder")
  {
    initialize();
  }

  Cylindrical::Cylindrical(const std::string& confstr) :
    BiorthBasis(confstr, "cylinder")
  {
    initialize();
  }


  double Cylindrical::DiskDens(double R, double z, double phi)
  {
    double ans = 0.0;

    switch (DTYPE) {

    case DiskType::constant:
      if (R < acyl && fabs(z) < hcyl)
	ans = 1.0/(2.0*hcyl*M_PI*acyl*acyl);
      break;
      
    case DiskType::gaussian:
      if (fabs(z) < hcyl)
	ans = 1.0/(2.0*hcyl*2.0*M_PI*acyl*acyl)*
	  exp(-R*R/(2.0*acyl*acyl));
      break;
      
    case DiskType::mn:
      {
	double Z2 = z*z + hcyl*hcyl;
	double Z  = sqrt(Z2);
	double Q2 = (acyl + Z)*(acyl + Z);
	ans = 0.25*hcyl*hcyl/M_PI*(acyl*R*R + (acyl + 3.0*Z)*Q2)/( pow(R*R + Q2, 2.5) * Z*Z2 );
      }
      break;
      
    case DiskType::doubleexpon:
      {
	double a1 = acyl;
	double a2 = acyl*aratio;
	double h1 = hcyl;
	double h2 = hcyl*hratio;
	double w1 = 1.0/(1.0+dweight);
	double w2 = dweight/(1.0+dweight);
	
	if (sech2) { h1 *= 0.5; h2 *= 0.5; }

	double f1 = cosh(z/h1);
	double f2 = cosh(z/h2);
	
	ans =
	  w1*exp(-R/a1)/(4.0*M_PI*a1*a1*h1*f1*f1) +
	  w2*exp(-R/a2)/(4.0*M_PI*a2*a2*h2*f2*f2) ;
      }
      break;

    case DiskType::diskbulge:
      {
	double h  = sech2 ? 0.5*hcyl : hcyl;
	double f  = cosh(z/h);
	double rr = pow(pow(R, 2) + pow(z,2), 0.5);
	double w1 = Mfac;
	double w2 = 1.0 - Mfac;
	double as = HERNA;
	
	ans = w1*exp(-R/acyl)/(4.0*M_PI*acyl*acyl*h*f*f) + 
	  w2*pow(as, 4)/(2.0*M_PI*rr)*pow(rr+as,-3.0) ;
      }
      break;
    case DiskType::python:
      ans = (*pyDens)(R, z, phi);
      break;
    case DiskType::exponential:
    default:
      {
	double h = sech2 ? 0.5*hcyl : hcyl;
	double f = cosh(z/h);
	ans = exp(-R/acyl)/(4.0*M_PI*acyl*acyl*h*f*f);
      }
      break;
    }
    
    if (rwidth>0.0) ans *= erf((rtrunc-R)/rwidth);
    
    return ans;
  }

  double Cylindrical::dcond(double R, double z, double phi, int M)
  {
    //
    // No shift for M==0
    //
    if (M==0) return DiskDens(R, z, phi);
    
    //
    // Fold into [-PI/M, PI/M] for M>=1
    //
    double dmult = M_PI/M, phiS;
    if (phi>M_PI)
      phiS = phi + dmult*(int)((2.0*M_PI - phi)/dmult);
    else
      phiS = phi - dmult*(int)(phi/dmult);
    
    //
    // Apply a shift along the x-axis
    //
    double x = R*cos(phiS) - ashift*acyl;
    double y = R*sin(phiS);
    
    return DiskDens(sqrt(x*x + y*y), z, atan2(y, x));
  }


  void Cylindrical::initialize()
  {
    // Basis identifier
    //
    BasisID = "Cylindrical";

    // Assign some defaults
    //
    rcylmin     = 0.001;
    rcylmax     = 20.0;
    acyl        = 0.01;
    hcyl        = 0.002;
    nmax        = 18;
    mmax        = 6;
    mlim        = std::numeric_limits<int>::max();
    lmaxfid     = 72;
    nmaxfid     = 64;
    ncylnx      = 256;
    ncylny      = 128;
    ncylodd     = 9;
    ncylr       = 2000;
    cachename   = ".eof_cache_file";
    oldcache    = false;
    Ignore      = false;
    deproject   = false;
    
    rnum        = 200;
    pnum        = 1;
    tnum        = 80;
    ashift      = 0.0;
    logarithmic = false;
    density     = true;
    EVEN_M      = false;
    cmapR       = 1;
    cmapZ       = 1;
    mtype       = "exponential";
    dtype       = "exponential";
    vflag       = 0;
    
    // Basis construction defaults
    //
    // The EmpCylSL deprojection from specified disk model (EXP or MN) -- change only if using MN
    dmodel      = "EXP"; 

    // Radial scale length ratio for disk basis construction with doubleexpon
    aratio      = 1.0; 

    // Vertical scale height ratio for disk basis construction with doubleexpon
    hratio      = 1.0;              

    // Ratio of second disk relative to the first disk for disk basis construction with double-exponential
    dweight     = 1.0;              

    // mass fraction for disk for diskbulge
    Mfac      = 1.0; 

    // Hernquist scale a disk basis construction with diskbulge
    HERNA      = 0.10; 

    // Width for erf truncation for EOF conditioning density (ignored if zero)
    rwidth      = 0.0;             

    // Fraction of scale length for shift in conditioning function
    ashift      = 0.0;              

    // Disk radial scaling factor for spherical deprojection model
    rfactor     = 1.0;  

    // Maximum disk radius for erf truncation of EOF conditioning density
    rtrunc      = 0.1;           

    // Power-law index for power-law disk profile
    ppow        = 4.0;

    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::Cylindrical", "parameter", unmatched, __FILE__, __LINE__);
    
    // Assign values from YAML
    //
    try {
      // These first two should not be user settable . . . but need them for now
      //
      if (conf["rcylmin"   ])    rcylmin  = conf["rcylmin"   ].as<double>();
      if (conf["rcylmax"   ])    rcylmax  = conf["rcylmax"   ].as<double>();
      
      if (conf["acyl"      ])       acyl  = conf["acyl"      ].as<double>();
      if (conf["hcyl"      ])       hcyl  = conf["hcyl"      ].as<double>();
      if (conf["sech2"     ])      sech2  = conf["sech2"     ].as<bool>();
      if (conf["lmaxfid"   ])    lmaxfid  = conf["lmaxfid"   ].as<int>();
      if (conf["nmaxfid"   ])    nmaxfid  = conf["nmaxfid"   ].as<int>();
      if (conf["nmax"      ])       nmax  = conf["nmax"      ].as<int>();
      if (conf["mmax"      ])       mmax  = conf["mmax"      ].as<int>();
      if (conf["mlim"      ])       mlim  = conf["mlim"      ].as<int>();
      if (conf["ncylnx"    ])     ncylnx  = conf["ncylnx"    ].as<int>();
      if (conf["ncylny"    ])     ncylny  = conf["ncylny"    ].as<int>();
      if (conf["ncylr"     ])      ncylr  = conf["ncylr"     ].as<int>();
      if (conf["ncylodd"   ])    ncylodd  = conf["ncylodd"   ].as<int>();
      if (conf["cachename" ])  cachename  = conf["cachename" ].as<std::string>();
      if (conf["eof_file"  ])  cachename  = conf["eof_file"  ].as<std::string>();
      if (conf["oldcache"  ])   oldcache  = conf["oldcache"  ].as<bool>();
      if (conf["rnum"      ])       rnum  = conf["rnum"      ].as<int>();
      if (conf["pnum"      ])       pnum  = conf["pnum"      ].as<int>();
      if (conf["tnum"      ])       tnum  = conf["tnum"      ].as<int>();

      if (conf["ashift"    ])     ashift  = conf["ashift"    ].as<double>();
      if (conf["logr"      ]) logarithmic = conf["logr"      ].as<bool>();
      if (conf["EVEN_M"    ])     EVEN_M  = conf["EVEN_M"    ].as<bool>();
      if (conf["cmapr"     ])      cmapR  = conf["cmapr"     ].as<int>();
      if (conf["cmapz"     ])      cmapZ  = conf["cmapz"     ].as<int>();
      if (conf["ignore"    ])      Ignore = conf["ignore"    ].as<bool>();
      if (conf["deproject" ])   deproject = conf["deproject" ].as<bool>();
      if (conf["dmodel"    ])      dmodel = conf["dmodel"    ].as<bool>();

      if (conf["aratio"    ])      aratio = conf["aratio"    ].as<double>();
      if (conf["hratio"    ])      hratio = conf["hratio"    ].as<double>();
      if (conf["dweight"   ])      dweight = conf["dweight"  ].as<double>();
      if (conf["Mfac"      ])      Mfac   = conf["Mfac"      ].as<double>();
      if (conf["HERNA"     ])      HERNA  = conf["HERNA"     ].as<double>();
      if (conf["rwidth"    ])      rwidth = conf["rwidth"    ].as<double>();
      if (conf["ashift"    ])      ashift = conf["ashift"    ].as<double>();
      if (conf["rfactor"   ])     rfactor = conf["rfactor"   ].as<double>();
      if (conf["rtrunc"    ])      rtrunc = conf["rtrunc"    ].as<double>();
      if (conf["pow"       ])        ppow = conf["ppow"      ].as<double>();
      if (conf["mtype"     ])       mtype = conf["mtype"     ].as<std::string>();
      if (conf["dtype"     ])       dtype = conf["dtype"     ].as<std::string>();
      if (conf["vflag"     ])       vflag = conf["vflag"     ].as<int>();
      if (conf["pyname"    ])      pyname = conf["pyname"    ].as<std::string>();
      if (conf["pcavar"]    )      pcavar = conf["pcavar"    ].as<bool>();

      // Deprecation warning
      if (conf["density"   ]) {
	if (myid==0)
	  std::cout << "Cylindrical: parameter 'density' is deprecated. "
		    << "The density field will be computed regardless."
		    << std::endl;
      }

      // Deprecation warning
      if (conf["eof_file"]) {
	if (myid==0)
	  std::cout << "Cylinder: parameter 'eof_file' is deprecated. "
		    << "and will be removed in a future release. Please "
		    << "use 'cachename' instead."
		    << std::endl;

	conf["cachename"] = conf["eof_file"];
      }
      
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing 'force' for Component <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Cylindrical: error parsing YAML");
    }
    
    // Enforce sane values for EOF integration
    //
    rnum = std::max<int>(10, rnum);
    pnum = std::max<int>(1,  pnum);
    tnum = std::max<int>(10, tnum);
    
    EmpCylSL::RMIN        = rcylmin;
    EmpCylSL::RMAX        = rcylmax;
    EmpCylSL::NUMX        = ncylnx;
    EmpCylSL::NUMY        = ncylny;
    EmpCylSL::NUMR        = ncylr;
    EmpCylSL::CMAPR       = cmapR;
    EmpCylSL::CMAPZ       = cmapZ;
    EmpCylSL::logarithmic = logarithmic;
    EmpCylSL::VFLAG       = vflag;
    
    // Check for non-null cache file name.  This must be specified
    // to prevent recomputation and unexpected behavior.
    //
    if (not conf["cachename"]) {
      throw std::runtime_error
	("Cylindrical requires a specified 'cachename' in your YAML config\n"
	 "for consistency with previous invocations and existing coefficient\n"
	 "sets.  Please add explicitly add 'cachename: name' to your config\n"
	 "with new 'name' for creating a basis or an existing 'name' for\n"
	 "reading a previously generated basis cache\n");
    }

    // Make the empirical orthogonal basis instance
    //
    sl = std::make_shared<EmpCylSL>
      (nmaxfid, lmaxfid, mmax, nmax, acyl, hcyl, ncylodd, cachename);
    
    // Set azimuthal harmonic order restriction?
    //
    if (mlim>=0)  sl->set_mlim(mlim);
    if (EVEN_M)   sl->setEven(EVEN_M);
      
    // Cache override for old Eigen cache
    //
    if (oldcache) sl->AllowOldCache();
    

    // Attempt to read EOF cache
    //
    if (sl->read_cache() == 0) {

      // Remake cylindrical basis
      //

      if (Ignore and myid==0) {
	std::cout << "---- BasisFactory: We have deprecated the 'ignore' parameter for the" << std::endl
		  << "----               Cylindrical class. If the cache file exists but does" << std::endl
		  << "----               not match the requested parameters, the old cache file" << std::endl
		  << "----               will be moved to a .bak file, a new basis will be re-" << std::endl
		  << "----               computed, and a new cache saved in its place.  Please" << std::endl
		  << "----               remove 'ignore' from your YAML configuration." << std::endl;
      }

      // Convert mtype string to lower case
      //
      std::transform(mtype.begin(), mtype.end(), mtype.begin(),
		     [](unsigned char c){ return std::tolower(c); });
      
      // Set EmpCylSL mtype.  This is the spherical function used to
      // generate the EOF basis.  If "deproject" is set, this will be
      // overriden in EmpCylSL.
      //
      EmpCylSL::mtype = EmpCylSL::Exponential;
      if (mtype.compare("exponential")==0)
	EmpCylSL::mtype = EmpCylSL::Exponential;
      else if (mtype.compare("gaussian")==0)
	EmpCylSL::mtype = EmpCylSL::Gaussian;
      else if (mtype.compare("plummer")==0)
	EmpCylSL::mtype = EmpCylSL::Plummer;
      else if (mtype.compare("power")==0) {
	EmpCylSL::mtype = EmpCylSL::Power;
	EmpCylSL::PPOW  = ppow;
      } else {
	if (myid==0) std::cout << "No EmpCylSL EmpModel named <"
			       << mtype << ">, valid types are: "
			       << "Exponential, Gaussian, Plummer, Power "
			       << "(not case sensitive)" << std::endl;
	throw std::runtime_error("Cylindrical:initialize: EmpCylSL bad parameter");
      }

      // Convert dtype string to lower case
      //
      std::transform(dtype.begin(), dtype.end(), dtype.begin(),
		     [](unsigned char c){ return std::tolower(c); });

      // Set DiskType.  This is the functional form for the disk used to
      // condition the basis.
      //
      try {				// Check for map entry, will through if the 
	DTYPE = dtlookup.at(dtype);	// key is not in the map.
	
	if (myid==0) {		// Report DiskType
	  std::cout << "---- DiskType is <" << dtype << ">" << std::endl;

	  if (not sech2) {
	    switch (DTYPE) {
	    case DiskType::doubleexpon:
	    case DiskType::exponential:
	    case DiskType::diskbulge:
	      std::cout << "---- pyEXP uses sech^2(z/h) rather than the more common sech^2(z/(2h))" << std::endl
			<< "---- Use the 'sech2: true' in your YAML config to use sech^2(z/(2h))" << std::endl
			<< "---- pyEXP will assume sech^2(z/(2h)) by default in v 7.9.0 and later" << std::endl;
	    }
	  }
	}
      }
      catch (const std::out_of_range& err) {
	if (myid==0) {
	  std::cout << "DiskType error in configuraton file" << std::endl;
	  std::cout << "Valid options are: ";
	  for (auto v : dtlookup) std::cout << v.first << " ";
	  std::cout << std::endl;
	}
	throw std::runtime_error("Cylindrical::initialize: invalid DiskType");
      }

      // Check for and initialize the Python density type
      //
      if (DTYPE == DiskType::python) {
	pyDens = std::make_shared<DiskDensityFunc>(pyname);
      }

      // Use these user models to deproject for the EOF spherical basis
      //
      if (deproject) {
	// The scale in EmpCylSL is assumed to be 1 so we compute the
	// height relative to the length
	//
	double H = sech2 ? 0.5*hcyl/acyl : hcyl/acyl;

	// The model instance (you can add others in DiskModels.H).
	// It's MN or Exponential if not MN.
	//
	EmpCylSL::AxiDiskPtr model;
	
	if (dmodel.compare("MN")==0) // Miyamoto-Nagai
	  model = std::make_shared<MNdisk>(1.0, H);
	else			// Default to exponential
	  model = std::make_shared<Exponential>(1.0, H);

	if (rwidth>0.0) {
	  model = std::make_shared<Truncated>(rtrunc/acyl,
					      rwidth/acyl,
					      model);
	  if (myid==0)
	    std::cout << "Made truncated model with R=" << rtrunc/acyl
		      << " and W=" << rwidth/acyl << std::endl;
	}
     
	sl->create_deprojection(H, rfactor, rnum, ncylr, model);
      }
    
      // Regenerate EOF from analytic density
      //
      std::function<double(double, double, double, int)> 
	f = [this](double R, double z, double phi, int M) -> double
	{
	  return this->dcond(R, z, phi, M);
	};

      sl->generate_eof(rnum, pnum, tnum, f);
    }

    // Orthogonality sanity check
    //
    if (myid==0) orthoTest();

    // Set cylindrical coordindates
    //
    coordinates = Coord::Cylindrical;
  }

  
  // Evaluate in spherical coordinates
  std::vector<double> Cylindrical::sph_eval(double r, double cth, double phi)
  {
    double sth = sqrt(1.0 - cth*cth);
    double R   = r*sth;
    double z   = r*cth;
    
    double tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp;
    
    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);

    double tpotr = tpotR*R/r + tpotz*z/R ;
    double tpott = tpotR*z/r - tpotz*R/r ;

    return
      {tdens0, tdens-tdens0, tdens,
       tpotl0, tpotl-tpotl0, tpotl, tpotr, tpott, tpotp};
  }
  
  // Evaluate in cartesian coordinates
  std::vector<double> Cylindrical::crt_eval(double x, double y, double z)
  {
    double R = sqrt(x*x + y*y);
    double phi = atan2(y, x);

    double tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp;

    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);

    double tpotx = tpotR*x/R - tpotp*y/R ;
    double tpoty = tpotR*y/R + tpotp*x/R ;

    return
      {tdens0, tdens - tdens0, tdens,
       tpotl0, tpotl - tpotl0, tpotl, tpotx, tpoty, tpotz};
  }
  
  // Evaluate in cartesian coordinates
  void Cylindrical::computeAccel(double x, double y, double z,
				 Eigen::Ref<Eigen::Vector3d> acc)
  {
    double R = sqrt(x*x + y*y);
    double phi = atan2(y, x);

    double tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp;

    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);

    double tpotx = tpotR*x/R - tpotp*y/R ;
    double tpoty = tpotR*y/R + tpotp*x/R ;

    acc << tpotx, tpoty, tpotz;
  }

  // Evaluate in cylindrical coordinates
  std::vector<double> Cylindrical::cyl_eval(double R, double z, double phi)
  {
    double tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp, height;

    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);

    if (midplane) {
      height = sl->accumulated_midplane_eval(R, -colh*hcyl, colh*hcyl, phi);
      return
	{tdens0, tdens - tdens0, tdens,
	 tpotl0, tpotl - tpotl0, tpotl,
	 tpotR, tpotz, tpotp, height};
    } else {
      return
	{tdens0, tdens - tdens0, tdens,
	 tpotl0, tpotl - tpotl0, tpotl,
	 tpotR, tpotz, tpotp};
    }
  }
  
  void Cylindrical::accumulate(double x, double y, double z, double mass)
  {
    double R   = sqrt(x*x + y*y);
    double phi = atan2(y, x);
    sl->accumulate(R, z, phi, mass, 0, 0);
  }
  
  void Cylindrical::reset_coefs(void)
  {
    sl->setup_accumulation();
  }
  
  void Cylindrical::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());

    cf->mmax   = mmax;
    cf->nmax   = nmax;
    cf->time   = time;

    Eigen::VectorXd cos1(nmax), sin1(nmax);

    // Initialize the values
    cos1.setZero();
    sin1.setZero();

    // Allocate the coefficient storage
    cf->store.resize((mmax+1)*nmax);

    // Create a new instance
    cf->coefs = std::make_shared<CoefClasses::CylStruct::coefType>
      (cf->store.data(), mmax+1, nmax);

    for (int m=0; m<=mmax; m++) {
      sl->get_coefs(m, cos1, sin1);

      for (int n=0; n<nmax; n++) {
	(*cf->coefs)(m, n) = {cos1(n), sin1(n)};
      }
    }
  }

  void Cylindrical::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    if (not dynamic_cast<CoefClasses::CylStruct*>(coef.get()))
      throw std::runtime_error("Cylindrical::set_coefs: you must pass a CoefClasses::CylStruct");

    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());

    // Cache the current coefficient structure
    //
    coefret = coef;

    for (int m=0; m<=mmax; m++) { // Set to zero on m=0 call only--------+
      sl->set_coefs(m, (*cf->coefs).row(m).real(), (*cf->coefs).row(m).imag(), m==0);
    }

    // Assign center if need be
    //
    if (cf->ctr.size())
      coefctr = cf->ctr;
    else
      coefctr = {0.0, 0.0, 0.0};
  }

  void Cylindrical::make_coefs(void)
  {
    sl->make_coefficients();
  }
  
  
  Cylindrical::BasisArray Cylindrical::getBasis
  (double xmin, double xmax, int numR, double zmin, double zmax, int numZ,
   bool linear)
  {
    // Allocate storage
    BasisArray ret(mmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numR, numZ); // Potential
	u["density"  ].resize(numR, numZ); // Density
	u["rforce"   ].resize(numR, numZ); // Radial force
	u["zforce"   ].resize(numR, numZ); // Vertical force
      }
    }
    
    // Grid spacing
    double delR = (xmax - xmin)/std::max<int>(numR-1, 1);
    double delZ = (zmax - zmin)/std::max<int>(numZ-1, 1);

    // Return values
    double p, d, fr, fz, fp;

    // Now, evaluate the grid
    for (int m=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	for (int i=0; i<numR; i++) {
	  double R = xmin + delR*i;
	  if (not linear) R = pow(10.0, R);
	  for (int j=0; j<numZ; j++) {
	    double Z = zmin + delZ*j;
	    sl->get_all(m, n, R, Z, 0.0, p, d, fr, fz, fp);
	    ret[m][n]["potential"](i,j) = p;
	    ret[m][n]["density"  ](i,j) = d;
	    ret[m][n]["rforce"   ](i,j) = fr;
	    ret[m][n]["zforce"   ](i,j) = fz;
	  }
	}
      }
    }

    return ret;
  }

  const std::set<std::string>
  FlatDisk::valid_keys = {
    "nmaxfid",
    "rcylmin",
    "rcylmax",
    "numx",
    "numy",
    "numr",
    "NQDHT",
    "knots",
    "logr",
    "model",
    "biorth",
    "scale",
    "rmin",
    "rmax",
    "self_consistent",
    "NO_M0",
    "NO_M1",
    "EVEN_M",
    "M0_BACK",
    "M0_ONLY",
    "NO_MONO",
    "diskconf",
    "background",
    "ssfrac",
    "playback",
    "coefMaster",
    "Lmax",
    "Mmax",
    "nmax",
    "mmax",
    "mlim",
    "dof",
    "subsamp",
    "samplesz",
    "vtkfreq",
    "tksmooth",
    "tkcum",
    "tk_type",
    "cachename"
  };

  FlatDisk::FlatDisk(const YAML::Node& CONF) :
    BiorthBasis(CONF, "flatdisk")
  {
    initialize();
  }

  FlatDisk::FlatDisk(const std::string& confstr) :
    BiorthBasis(confstr, "flatdisk")
  {
    initialize();
  }

  void FlatDisk::initialize()
  {
    // Basis identifier
    //
    BasisID = "FlatDisk";

    // Assign some defaults
    //
    cmap       = 1;
    mmax       = 6;
    nmax       = 18;
    
    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::FlatDisk", "parameter", unmatched, __FILE__, __LINE__);
    
    // Default cachename, empty by default
    //
    std::string cachename;

    // Assign values from YAML
    //
    try {
      if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
      if (conf["Lmax"])      mmax       = conf["Lmax"].as<int>(); // Proxy
      if (conf["Mmax"])      mmax       = conf["Mmax"].as<int>(); // Proxy
      if (conf["mmax"])      mmax       = conf["mmax"].as<int>();
      if (conf["nmax"])      nmax       = conf["nmax"].as<int>();
      
      if (conf["rcylmin"]) 
	rcylmin = conf["rcylmin"].as<double>();
      else
	rcylmin = 0.0;
      
      if (conf["rcylmax"]) 
	rcylmax = conf["rcylmax"].as<double>();
      else
	rcylmax = 10.0;
      
      N1 = 0;
      N2 = std::numeric_limits<int>::max();
      NO_M0 = NO_M1 = EVEN_M = M0_only = false;
      
      if (conf["N1"]   )     N1        = conf["N1"].as<bool>();
      if (conf["N2"]   )     N2        = conf["N2"].as<bool>();
      if (conf["NO_M0"])     NO_M0     = conf["NO_M0"].as<bool>();
      if (conf["NO_M1"])     NO_M1     = conf["NO_M1"].as<bool>();
      if (conf["EVEN_M"])    EVEN_M    = conf["EVEN_M"].as<bool>();
      if (conf["M0_ONLY"])   M0_only   = conf["M0_ONLY"].as<bool>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("FlatDisk: error parsing YAML");
    }
    
    // Set cmapR and cmapZ defaults
    //
    if (not conf["cmapR"])   conf["cmapR"] = cmap;
    if (not conf["cmapZ"])   conf["cmapZ"] = cmap;
    
    // Set characteristic radius defaults
    //
    if (not conf["scale"])   conf["scale"]   = 1.0;

    // Check for non-null cache file name.  This must be specified
    // to prevent recomputation and unexpected behavior.
    //
    if (not conf["cachename"]) {
      throw std::runtime_error
	("FlatDisk requires a specified cachename in your YAML config\n"
	 "for consistency with previous invocations and existing coefficient\n"
	 "sets.  Please add explicitly add 'cachename: name' to your config\n"
	 "with new 'name' for creating a basis or an existing 'name' for\n"
	 "reading a previously generated basis cache\n");
    }

    // Finally, make the basis
    //
    ortho = std::make_shared<BiorthCyl>(conf);
    
    // Orthogonality sanity check
    //
    if (myid==0) orthoTest();

    // Get max threads
    //
    int nthrds = omp_get_max_threads();

    // Allocate memory
    //
    potd.resize(nthrds);
    potR.resize(nthrds);
    potZ.resize(nthrds);
    dend.resize(nthrds);

    for (auto & v : potd) v.resize(mmax+1, nmax);
    for (auto & v : potR) v.resize(mmax+1, nmax);
    for (auto & v : potZ) v.resize(mmax+1, nmax);
    for (auto & v : dend) v.resize(mmax+1, nmax);

    expcoef.resize(2*mmax+1, nmax);
    expcoef.setZero();
      
    work.resize(nmax);
      
    used = 0;

    // Set cylindrical coordindates
    //
    coordinates = Coord::Cylindrical;
  }
  
  void FlatDisk::reset_coefs(void)
  {
    if (expcoef.rows()>0 && expcoef.cols()>0) expcoef.setZero();
    totalMass = 0.0;
    used = 0;
  }
  
  
  void FlatDisk::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());

    cf->mmax   = mmax;
    cf->nmax   = nmax;
    cf->time   = time;

    // Allocate the coefficient storage
    cf->store.resize((mmax+1)*nmax);

    // Make the coefficient map
    cf->coefs = std::make_shared<CoefClasses::CylStruct::coefType>
      (cf->store.data(), mmax+1, nmax);

    for (int m=0, m0=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	if (m==0)
	  (*cf->coefs)(m, n) = {expcoef(m0, n), 0.0};
	else
	  (*cf->coefs)(m, n) = {expcoef(m0, n), expcoef(m0+1, n)};
      }
      if (m==0) m0 += 1;
      else      m0 += 2;
    }
  }

  void FlatDisk::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (not dynamic_cast<CoefClasses::CylStruct*>(coef.get()))
      throw std::runtime_error("FlatDisk::set_coefs: you must pass a CoefClasses::CylStruct");

    // Sanity check on dimensionality
    //
    {
      auto cc = dynamic_cast<CoefClasses::CylStruct*>(coef.get());
      auto cf = cc->coefs;
      int rows = cf->rows();
      int cols = cf->cols();
      if (rows != mmax+1 or cols != nmax) {
	std::ostringstream sout;
	sout << "FlatDisk::set_coefs: the basis has (mmax+1, nmax)=("
	     << mmax+1 << ", " << nmax
	     << "). The coef structure has (rows, cols)=("
	     << rows << ", " << cols << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());
    auto & cc = *cf->coefs;

    // Cache the current coefficient structure
    //
    coefret = coef;

    // Assign internal coefficient table (doubles) from the complex struct
    //
    for (int m=0, m0=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	if (m==0)
	  expcoef(m0,   n) = cc(m, n).real();
	else {
	  expcoef(m0,   n) = cc(m, n).real();
	  expcoef(m0+1, n) = cc(m, n).imag();
	}
      }
      if (m==0) m0 += 1;
      else      m0 += 2;
    }

    // Assign center if need be
    //
    if (cf->ctr.size())
      coefctr = cf->ctr;
    else
      coefctr = {0.0, 0.0, 0.0};
  }

  void FlatDisk::accumulate(double x, double y, double z, double mass)
  {
    // Normalization factors
    //
    constexpr double norm0 = 1.0;
    constexpr double norm1 = M_SQRT2;

    //======================
    // Compute coefficients 
    //======================
    
    double R2 = x*x + y*y;
    double R  = sqrt(R2);
    
    // Get thread id
    int tid = omp_get_thread_num();

    if (R < ortho->getRtable() and fabs(z) < ortho->getRtable()) {
    
      used++;
      totalMass += mass;
    
      double phi = atan2(y, x);

      ortho->get_pot(potd[tid], R, 0.0);
    
      // M loop
      for (int m=0, moffset=0; m<=mmax; m++) {
	
	if (m==0) {
	  for (int n=0; n<nmax; n++) {
	    expcoef(moffset, n) += potd[tid](m, n)* mass * norm0;
	  }
	  
	  moffset++;
	}
	else {
	  double ccos = cos(phi*m);
	  double ssin = sin(phi*m);
	  for (int n=0; n<nmax; n++) {
	    expcoef(moffset  , n) += ccos * potd[tid](m, n) * mass * norm1;
	    expcoef(moffset+1, n) += ssin * potd[tid](m, n) * mass * norm1;
	  }
	  moffset+=2;
	}
      }
    }
    
  }
  
  void FlatDisk::make_coefs()
  {
    if (use_mpi) {
      
      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      for (int m=0; m<2*mmax+1; m++) {
	work = expcoef.row(m);
	MPI_Allreduce(MPI_IN_PLACE, work.data(), nmax, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	expcoef.row(m) = work;
      }
    }
  }
  
  std::vector<double>FlatDisk::cyl_eval(double R, double z, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Fixed values
    constexpr double norm0 = 1.0;
    constexpr double norm1 = M_SQRT2;

    double den0=0, den1=0, pot0=0, pot1=0, rpot=0, zpot=0, ppot=0;

    // Off grid evaluation
    if (R>ortho->getRtable() or fabs(z)>ortho->getRtable()) {
      double r2 = R*R + z*z;
      double r  = sqrt(r2);
      pot0 = -totalMass/r;
      rpot = -totalMass*R/(r*r2 + 10.0*std::numeric_limits<double>::min());
      zpot = -totalMass*z/(r*r2 + 10.0*std::numeric_limits<double>::min());
      
      return {den0, den1, den0+den1, pot0, pot1, pot0+pot1, rpot, zpot, ppot};
    }

    // Get the basis fields
    //
    ortho->get_dens   (dend[tid],  R, z);
    ortho->get_pot    (potd[tid],  R, z);
    ortho->get_rforce (potR[tid],  R, z);
    ortho->get_zforce (potZ[tid],  R, z);
    
    // m loop
    //
    for (int m=0, moffset=0; m<=mmax; m++) {
      
      if (m==0 and NO_M0)        { moffset++;    continue; }
      if (m==1 and NO_M1)        { moffset += 2; continue; }
      if (EVEN_M and m/2*2 != m) { moffset += 2; continue; }
      if (m>0 and M0_only)       break;

      if (m==0) {
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  den0 += expcoef(0, n) * dend[tid](0, n) * norm0;
	  pot0 += expcoef(0, n) * potd[tid](0, n) * norm0;
	  rpot += expcoef(0, n) * potR[tid](0, n) * norm0;
	  zpot += expcoef(0, n) * potZ[tid](0, n) * norm0;
	}
	
	moffset++;
      } else {
	double cosm = cos(phi*m), sinm = sin(phi*m);
	double vc, vs;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * dend[tid](m, n);
	  vs += expcoef(moffset+1, n) * dend[tid](m, n);
	}
	
	den1 += (vc*cosm + vs*sinm) * norm1;
      
	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potd[tid](m, n);
	  vs += expcoef(moffset+1, n) * potd[tid](m, n);
	}
	
	pot1 += ( vc*cosm + vs*sinm) * norm1;
	ppot += (-vc*sinm + vs*cosm) * m * norm1;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potR[tid](m, n);
	  vs += expcoef(moffset+1, n) * potR[tid](m, n);
	}

	rpot += (vc*cosm + vs*sinm) * norm1;
	
	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potZ[tid](m, n);
	  vs += expcoef(moffset+1, n) * potZ[tid](m, n);
	}

	zpot += (vc*cosm + vs*sinm) * norm1;

	moffset +=2;
      }
    }

    den0 *= -1.0;
    den1 *= -1.0;
    pot0 *= -1.0;
    pot1 *= -1.0;
    rpot *= -1.0;
    zpot *= -1.0;
    ppot *= -1.0;

    return {den0, den1, den0+den1, pot0, pot1, pot0+pot1, rpot, zpot, ppot};
  }

  void FlatDisk::computeAccel(double x, double y, double z, 
			      Eigen::Ref<Eigen::Vector3d> acc)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Fixed values
    constexpr double norm0 = 0.5*M_2_SQRTPI/M_SQRT2;
    constexpr double norm1 = 0.5*M_2_SQRTPI;

    // Compute polar coordinates
    double R   = std::sqrt(x*x + y*y);
    double phi = std::atan2(y, x);
    

    double rpot=0, zpot=0, ppot=0;

    // Off grid evaluation
    if (R>ortho->getRtable() or fabs(z)>ortho->getRtable()) {
      double r2 = R*R + z*z;
      double r  = sqrt(r2);

      rpot = -totalMass*R/(r*r2 + 10.0*std::numeric_limits<double>::min());
      zpot = -totalMass*z/(r*r2 + 10.0*std::numeric_limits<double>::min());
      
      acc << rpot, zpot, ppot;
    }

    // Get the basis fields
    //
    ortho->get_pot    (potd[tid],  R, z);
    ortho->get_rforce (potR[tid],  R, z);
    ortho->get_zforce (potZ[tid],  R, z);
    
    // m loop
    //
    for (int m=0, moffset=0; m<=mmax; m++) {
      
      if (m==0 and NO_M0)        { moffset++;    continue; }
      if (m==1 and NO_M1)        { moffset += 2; continue; }
      if (EVEN_M and m/2*2 != m) { moffset += 2; continue; }
      if (m>0 and M0_only)       break;

      if (m==0) {
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  rpot += expcoef(0, n) * potR[tid](0, n) * norm0;
	  zpot += expcoef(0, n) * potZ[tid](0, n) * norm0;
	}
	
	moffset++;
      } else {
	double cosm = cos(phi*m), sinm = sin(phi*m);
	double vc, vs;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potd[tid](m, n);
	  vs += expcoef(moffset+1, n) * potd[tid](m, n);
	}
	
	ppot += (-vc*sinm + vs*cosm) * m * norm1;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potR[tid](m, n);
	  vs += expcoef(moffset+1, n) * potR[tid](m, n);
	}

	rpot += (vc*cosm + vs*sinm) * norm1;
	
	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potZ[tid](m, n);
	  vs += expcoef(moffset+1, n) * potZ[tid](m, n);
	}

	zpot += (vc*cosm + vs*sinm) * norm1;

	moffset +=2;
      }
    }

    rpot *= -1.0;
    zpot *= -1.0;
    ppot *= -1.0;

    double potx = rpot*x/R - ppot*y/R;
    double poty = rpot*y/R + ppot*x/R;

    acc << potx, poty, zpot;
  }


  std::vector<double> FlatDisk::sph_eval(double r, double costh, double phi)
  {
    // Cylindrical coords
    //
    double sinth = sqrt(fabs(1.0 - costh*costh));
    double R = r*sinth, z = r*costh;

    auto v = cyl_eval(R, z, phi);
  
    // Spherical force element converstion
    //
    double potr = v[6]*sinth + v[7]*costh;
    double pott = v[6]*costh - v[7]*sinth;

    return {v[0], v[1], v[2], v[3], v[4], v[5], potr, pott, v[8]};
  }

  std::vector<double> FlatDisk::crt_eval(double x, double y, double z)
  {
    // Cylindrical coords from Cartesian
    //
    double R = sqrt(x*x + y*y) + 1.0e-18;
    double phi = atan2(y, x);

    auto v = cyl_eval(R, z, phi);

    double potx = v[6]*x/R - v[8]*y/R;
    double poty = v[6]*y/R + v[8]*x/R;

    return {v[0], v[1], v[2], v[3], v[4], v[5], potx, poty, v[7]};
  }

  std::vector<Eigen::MatrixXd> FlatDisk::orthoCheck(int knots)
  {
    return ortho->orthoCheck();
  }
  

  FlatDisk::BasisArray FlatDisk::getBasis
  (double logxmin, double logxmax, int numgrid)
  {
    // Assing return storage
    BasisArray ret(mmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numgrid); // Potential
	u["density"  ].resize(numgrid); // Density
	u["rforce"   ].resize(numgrid); // Radial force
      }
    }

    // Radial grid spacing
    double dx = (logxmax - logxmin)/(numgrid-1);

    // Basis storage
    Eigen::MatrixXd tabpot, tabden, tabrfc;

    // Evaluate on the plane
    for (int i=0; i<numgrid; i++) {
      ortho->get_pot   (tabpot, pow(10.0, logxmin + dx*i), 0.0);
      ortho->get_dens  (tabden, pow(10.0, logxmin + dx*i), 0.0);
      ortho->get_rforce(tabrfc, pow(10.0, logxmin + dx*i), 0.0);
      for (int m=0; m<=mmax; m++) {
	for (int n=0; n<nmax; n++){
	  ret[m][n]["potential"](i) = tabpot(m, n);
	  ret[m][n]["density"  ](i) = tabden(m, n);
	  ret[m][n]["rforce"   ](i) = tabrfc(m, n);
	}
      }
    }
    
    return ret;
  }

  const std::set<std::string>
  CBDisk::valid_keys = {
    "self_consistent",
    "NO_M0",
    "NO_M1",
    "EVEN_M",
    "M0_BACK",
    "M0_ONLY",
    "NO_MONO",
    "background",
    "playback",
    "coefMaster",
    "scale",
    "Lmax",
    "Mmax",
    "nmax",
    "mmax",
    "mlim",
    "dof",
    "subsamp",
    "samplesz",
    "vtkfreq",
    "tksmooth",
    "tkcum",
    "tk_type"
  };

  CBDisk::CBDisk(const YAML::Node& CONF) :
    BiorthBasis(CONF, "CBDisk")
  {
    initialize();
  }

  CBDisk::CBDisk(const std::string& confstr) :
    BiorthBasis(confstr, "CBDisk")
  {
    initialize();
  }

  void CBDisk::initialize()
  {
    // Basis identifier
    //
    BasisID = "CBDisk";

    // Assign some defaults
    //
    mmax       = 6;
    nmax       = 18;
    scale      = 1.0;
    
    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::CBDisk", "parameter", unmatched, __FILE__, __LINE__);
    
    // Assign values from YAML
    //
    try {
      if (conf["Lmax"])      mmax       = conf["Lmax"].as<int>(); // Proxy
      if (conf["Mmax"])      mmax       = conf["Mmax"].as<int>(); // Proxy
      if (conf["mmax"])      mmax       = conf["mmax"].as<int>();
      if (conf["nmax"])      nmax       = conf["nmax"].as<int>();
      if (conf["scale"])     scale      = conf["scale"].as<double>();
      
      N1 = 0;
      N2 = std::numeric_limits<int>::max();
      NO_M0 = NO_M1 = EVEN_M = M0_only = false;
      
      if (conf["N1"]   )     N1        = conf["N1"].as<bool>();
      if (conf["N2"]   )     N2        = conf["N2"].as<bool>();
      if (conf["NO_M0"])     NO_M0     = conf["NO_M0"].as<bool>();
      if (conf["NO_M1"])     NO_M1     = conf["NO_M1"].as<bool>();
      if (conf["EVEN_M"])    EVEN_M    = conf["EVEN_M"].as<bool>();
      if (conf["M0_ONLY"])   M0_only   = conf["M0_ONLY"].as<bool>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("CBDisk: error parsing YAML");
    }
    
    // Set characteristic radius defaults
    //
    if (not conf["scale"])   conf["scale"]   = 1.0;

    // Potential, force, and density scaling
    //
    fac1 = pow(scale, -0.5);
    fac2 = pow(scale, -1.5);

    // Orthogonality sanity check
    //
    if (myid==0) orthoTest();

    // Get max threads
    //
    int nthrds = omp_get_max_threads();

    // Allocate memory
    //
    potd.resize(nthrds);
    potR.resize(nthrds);
    dend.resize(nthrds);

    for (auto & v : potd) v.resize(mmax+1, nmax);
    for (auto & v : potR) v.resize(mmax+1, nmax);
    for (auto & v : dend) v.resize(mmax+1, nmax);

    expcoef.resize(2*mmax+1, nmax);
    expcoef.setZero();
      
    work.resize(nmax);
      
    used = 0;

    // Set cylindrical coordindates
    //
    coordinates = Coord::Cylindrical;
  }
  
  // Get potential
  void CBDisk::get_pot(Eigen::MatrixXd& tab, double r)
  {
    tab.resize(mmax+1, nmax);
    Eigen::VectorXd a(nmax);
    for (int m=0; m<=mmax; m++) {
      potl(m, r/scale, a);
      tab.row(m) = a;
    }
    tab *= fac1;
  }

  // Get density
  void CBDisk::get_dens(Eigen::MatrixXd& tab, double r)
  {
    tab.resize(mmax+1, nmax);
    Eigen::VectorXd a(nmax);
    for (int m=0; m<=mmax; m++) {
      dens(m, r/scale, a);
      tab.row(m) = a;
    }
    tab *= fac2;
  }

  // Get force
  void CBDisk::get_force(Eigen::MatrixXd& tab, double r)
  {
    tab.resize(mmax+1, nmax);
    Eigen::VectorXd a(nmax);
    for (int m=0; m<=mmax; m++) {
      dpot(m, r/scale, a);
      tab.row(m) = a;
    }
    tab *= fac1/scale;
  }

  //  Routines for computing biorthonormal pairs based on 
  //  Clutton-Brock's 2-dimensional series
  //
  double CBDisk::phif(const int n, const int m, const double r)
  {
    // By recurrance relation
    //
    double r2 = r*r;
    double fac = 1.0/(1.0 + r2);
    double cur = sqrt(fac);

    for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);

    if (n==0) return cur;
  
    double curl1 = cur;
    double curl2;
  
    fac *= r2 - 1.0;
    cur *= fac*(2*m+1);

    if (n==1) return cur;

    for (int nn=2; nn<=n; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur   = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
	(1.0 + (double)(2*m-1)/nn)*curl2;
    }
    
    return cur;
  }

  double CBDisk::potl(const int n, const int m, const double r)
  {
    return pow(r, m) * phif(n, m, r)/sqrt(norm(n, m));
  }

  double CBDisk::dpot(const int n, const int m, const double r)
  {
    double ret = dphi(n, m, r);
    if (m) ret = (phif(n, m, r)*m/r + ret) * pow(r, m);
    return ret/sqrt(norm(n, m));
  }
  
  double CBDisk::dphi(const int n, const int m, const double r)
  {
    double ret = phif(n, m+1, r);
    if (n>0) ret -= 2.0*phif(n-1, m+1, r);
    if (n>1) ret += phif(n-2, m+1, r);
    return -r*ret;
  }


  // By recurrance relation
  //
  void CBDisk::potl(const int m, const double r, Eigen::VectorXd& a)
  {
    a.resize(nmax);

    double pfac = pow(r, m);
  
    double r2  = r*r;
    double fac = 1.0/(1.0 + r2);
    double cur = sqrt(fac);
    
    for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);
    
    a(0) = pfac*cur;
    
    if (nmax>0) {
      
      double curl1 = cur;
      double curl2;
      
      fac *= r2 - 1.0;
      cur *= fac*(2*m+1);
      
      a(1) = pfac*cur;
      
      for (int nn=2; nn<nmax; nn++) {
	curl2 = curl1;
	curl1 = cur;
	cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
	  (1.0 + (double)(2*m-1)/nn)*curl2;
	a(nn) = pfac*cur;
      }
    }
    
    for (int n=0; n<nmax; n++) a(n) /= sqrt(norm(n, m));
    
    return;
  }

  // By recurrance relation
  //
  void CBDisk::dpot(const int m, const double r, Eigen::VectorXd& a)
  {
    a.resize(nmax);
    potl(m, r, a);
    a *= static_cast<double>(m)/r;

    Eigen::VectorXd b(nmax);
    potl(m+1, r, b);

    for (int n=0; n<nmax; n++) {
      a(n) -= b(n);
      if (n>0) a(n) += 2.0*b(n-1);
      if (n>1) a(n) -= b(n-2);
    }
  }
  
  double CBDisk::dens(int n, int m, double r)
  {
    double f = 0.5/sqrt(norm(n, m))/M_PI;
    
    if (n>=2) 
      return f*pow(r, m)*(phif(n, m+1, r) - phif(n-2, m+1, r));
    else
      return f*pow(r, m)*phif(n, m+1, r);
  }

  
  void CBDisk::dens(int mm, double r, Eigen::VectorXd& a)
  {
    a.resize(nmax);
  
    double pfac = pow(r, (double)mm+1.0e-20);
    
    int m = mm + 1;
    double r2 = r*r;
    double fac = 1.0/(1.0 + r2);
    double cur = sqrt(fac);
    
    for (int M=1; M<=m; M++) cur *= fac*(2*M - 1);

    a(0) = pfac*cur;

    if (nmax>0) {
  
      double curl1 = cur;
      double curl2;
      
      fac *= r2 - 1.0;
      cur *= fac*(2*m+1);
      
      a(1) = pfac*cur;
      
      for (int nn=2; nn<nmax; nn++) {
	curl2 = curl1;
	curl1 = cur;
	cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
	  (1.0 + (double)(2*m-1)/nn)*curl2;
	a(nn) = pfac*cur;
      }
      
      for (int nn=nmax-1; nn>1; nn--)
	a(nn) -= a(nn-2);
    }
  
    for (int n=0; n<nmax; n++) a(n) *= 0.5/sqrt(norm(n, mm))/M_PI;

    return;
  }
  

  double CBDisk::norm(int n, int m)
  {
  double ans = 1.0;
  
  for (int i=n+1; i<=n+2*m; i++) ans *= i;

  return pow(0.5, 2*m+1)*ans;
  }
  
  std::vector<Eigen::MatrixXd> CBDisk::orthoCheck(int num)
  {
    const double tol = 1.0e-4;
    LegeQuad lq(num);

    std::vector<Eigen::MatrixXd> ret(mmax+1);
    for (auto & v : ret) v.resize(nmax, nmax);

    double Rmax = scale*100.0;

    for (int m=0; m<=mmax; m++) {

      ret[m].setZero();

      for (int i=0; i<num; i++) {
	double r = lq.knot(i) * Rmax, fac = lq.weight(i) * Rmax;
	Eigen::VectorXd vpot, vden;
      
	potl(m, r/scale, vpot);
	dens(m, r/scale, vden);

	for (int j=0; j<nmax; j++) {
	  for (int k=0; k<nmax; k++) {
	    ret[m](j, k) += fac * vpot(j) * vden(k) * r * 2.0*M_PI * fac1*fac2;
	  }
	}
      }
    }

    // DEBUG
    std::ofstream out("orthoCheck.dat");
    out << "scale=" << scale << " fac1=" << fac1 << " fac2=" << fac2 << std::endl;
    for (int m=0; m<=mmax; m++) {
      out << std::string(80, '-') << std::endl
	  << "---- m=" << m       << std::endl
	  << std::string(80, '-') << std::endl
	  << ret[m]               << std::endl
	  << std::string(80, '-') << std::endl << std::endl;
    }
    // END DEBUG

    return ret;
  }

  void CBDisk::reset_coefs(void)
  {
    if (expcoef.rows()>0 && expcoef.cols()>0) expcoef.setZero();
    totalMass = 0.0;
    used = 0;
  }
  
  
  void CBDisk::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());

    cf->mmax   = mmax;
    cf->nmax   = nmax;
    cf->time   = time;

    // Allocate the coefficient storage
    cf->store.resize((mmax+1)*nmax);

    // Make the coefficient map
    cf->coefs = std::make_shared<CoefClasses::CylStruct::coefType>
      (cf->store.data(), mmax+1, nmax);

    for (int m=0, m0=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	if (m==0)
	  (*cf->coefs)(m, n) = {expcoef(m0, n), 0.0};
	else
	  (*cf->coefs)(m, n) = {expcoef(m0, n), expcoef(m0+1, n)};
      }
      if (m==0) m0 += 1;
      else      m0 += 2;
    }
  }

  void CBDisk::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (not dynamic_cast<CoefClasses::CylStruct*>(coef.get()))
      throw std::runtime_error("CBDisk::set_coefs: you must pass a CoefClasses::CylStruct");

    // Sanity check on dimensionality
    //
    {
      auto cc = dynamic_cast<CoefClasses::CylStruct*>(coef.get());
      auto cf = cc->coefs;
      int rows = cf->rows();
      int cols = cf->cols();
      if (rows != mmax+1 or cols != nmax) {
	std::ostringstream sout;
	sout << "CBDisk::set_coefs: the basis has (mmax+1, nmax)=("
	     << mmax+1 << ", " << nmax
	     << "). The coef structure has (rows, cols)=("
	     << rows << ", " << cols << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());
    auto & cc = *cf->coefs;

    // Cache the current coefficient structure
    coefret = coef;

    // Assign internal coefficient table (doubles) from the complex struct
    //
    for (int m=0, m0=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	if (m==0)
	  expcoef(m0,   n) = cc(m, n).real();
	else {
	  expcoef(m0,   n) = cc(m, n).real();
	  expcoef(m0+1, n) = cc(m, n).imag();
	}
      }
      if (m==0) m0 += 1;
      else      m0 += 2;
    }

    // Assign center if need be
    //
    if (cf->ctr.size())
      coefctr = cf->ctr;
    else
      coefctr = {0.0, 0.0, 0.0};
  }

  void CBDisk::accumulate(double x, double y, double z, double mass)
  {
    // Normalization factors
    //
    constexpr double norm0 = 1.0;
    constexpr double norm1 = M_SQRT2;

    //======================
    // Compute coefficients 
    //======================
    
    double R2 = x*x + y*y;
    double R  = sqrt(R2);
    
    // Get thread id
    int tid = omp_get_thread_num();

    used++;
    totalMass += mass;
    
    double phi = atan2(y, x);

    get_pot(potd[tid], R);
    
    // M loop
    for (int m=0, moffset=0; m<=mmax; m++) {
	
      if (m==0) {
	for (int n=0; n<nmax; n++) {
	  expcoef(moffset, n) += potd[tid](m, n)* mass * norm0;
	}
	
	moffset++;
      }
      else {
	double ccos = cos(phi*m);
	double ssin = sin(phi*m);
	for (int n=0; n<nmax; n++) {
	  expcoef(moffset  , n) += ccos * potd[tid](m, n) * mass * norm1;
	  expcoef(moffset+1, n) += ssin * potd[tid](m, n) * mass * norm1;
	}
	moffset+=2;
      }
    }
  }
  
  void CBDisk::make_coefs()
  {
    if (use_mpi) {
      
      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      for (int m=0; m<2*mmax+1; m++) {
	work = expcoef.row(m);
	MPI_Allreduce(MPI_IN_PLACE, work.data(), nmax, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	expcoef.row(m) = work;
      }
    }
  }
  
  std::vector<double> CBDisk::cyl_eval(double R, double z, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Fixed values
    constexpr double norm0 = 1.0;
    constexpr double norm1 = M_SQRT2;

    double den0=0, den1=0, pot0=0, pot1=0, rpot=0, zpot=0, ppot=0;

    // Get the basis fields
    //
    get_dens  (dend[tid], R);
    get_pot   (potd[tid], R);
    get_force (potR[tid], R);
    
    // m loop
    //
    for (int m=0, moffset=0; m<=mmax; m++) {
      
      if (m==0 and NO_M0)        { moffset++;    continue; }
      if (m==1 and NO_M1)        { moffset += 2; continue; }
      if (EVEN_M and m/2*2 != m) { moffset += 2; continue; }
      if (m>0 and M0_only)       break;

      if (m==0) {
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  den0 += expcoef(0, n) * dend[tid](0, n) * norm0;
	  pot0 += expcoef(0, n) * potd[tid](0, n) * norm0;
	  rpot += expcoef(0, n) * potR[tid](0, n) * norm0;
	}
	
	moffset++;
      } else {
	double cosm = cos(phi*m), sinm = sin(phi*m);
	double vc, vs;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * dend[tid](m, n);
	  vs += expcoef(moffset+1, n) * dend[tid](m, n);
	}
	
	den1 += (vc*cosm + vs*sinm) * norm1;
      
	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potd[tid](m, n);
	  vs += expcoef(moffset+1, n) * potd[tid](m, n);
	}
	
	pot1 += ( vc*cosm + vs*sinm) * norm1;
	ppot += (-vc*sinm + vs*cosm) * m * norm1;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potR[tid](m, n);
	  vs += expcoef(moffset+1, n) * potR[tid](m, n);
	}

	rpot += (vc*cosm + vs*sinm) * norm1;
	
	moffset +=2;
      }
    }

    den0 *= -1.0;
    den1 *= -1.0;
    pot0 *= -1.0;
    pot1 *= -1.0;
    rpot *= -1.0;
    ppot *= -1.0;

    return {den0, den1, den0+den1, pot0, pot1, pot0+pot1, rpot, zpot, ppot};
  }


  void CBDisk::computeAccel(double x, double y, double z,
			    Eigen::Ref<Eigen::Vector3d> acc)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Fixed values
    constexpr double norm0 = 1.0;
    constexpr double norm1 = M_SQRT2;

    double R   = std::sqrt(x*x + y*y);
    double phi = std::atan2(y, x);

    double rpot=0, zpot=0, ppot=0;

    // Get the basis fields
    //
    get_pot   (potd[tid], R);
    get_force (potR[tid], R);
    
    // m loop
    //
    for (int m=0, moffset=0; m<=mmax; m++) {
      
      if (m==0 and NO_M0)        { moffset++;    continue; }
      if (m==1 and NO_M1)        { moffset += 2; continue; }
      if (EVEN_M and m/2*2 != m) { moffset += 2; continue; }
      if (m>0 and M0_only)       break;

      if (m==0) {
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  rpot += expcoef(0, n) * potR[tid](0, n) * norm0;
	}
	
	moffset++;
      } else {
	double cosm = cos(phi*m), sinm = sin(phi*m);
	double vc, vs;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potd[tid](m, n);
	  vs += expcoef(moffset+1, n) * potd[tid](m, n);
	}
	
	ppot += (-vc*sinm + vs*cosm) * m * norm1;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potR[tid](m, n);
	  vs += expcoef(moffset+1, n) * potR[tid](m, n);
	}

	rpot += (vc*cosm + vs*sinm) * norm1;
	
	moffset +=2;
      }
    }

    rpot *= -1.0;
    ppot *= -1.0;


    double potx = rpot*x/R - ppot*y/R;
    double poty = rpot*y/R + ppot*x/R;

    acc << potx, poty, zpot;
  }


  std::vector<double> CBDisk::sph_eval(double r, double costh, double phi)
  {
    // Cylindrical coords
    //
    double sinth = sqrt(fabs(1.0 - costh*costh));
    double R = r*sinth, z = r*costh;

    auto v = cyl_eval(R, z, phi);
  
    // Spherical force element converstion
    //
    double potr = v[6]*sinth + v[7]*costh;
    double pott = v[6]*costh - v[7]*sinth;

    return {v[0], v[1], v[2], v[3], v[4], v[5], potr, pott, v[8]};
  }

  std::vector<double> CBDisk::crt_eval(double x, double y, double z)
  {
    // Cylindrical coords from Cartesian
    //
    double R = sqrt(x*x + y*y) + 1.0e-18;
    double phi = atan2(y, x);

    auto v = cyl_eval(R, z, phi);

    double potx = v[6]*x/R - v[8]*y/R;
    double poty = v[6]*y/R + v[8]*x/R;

    return {v[0], v[1], v[2], v[3], v[4], v[5], potx, poty, v[7]};
  }

  CBDisk::BasisArray CBDisk::getBasis
  (double logxmin, double logxmax, int numgrid)
  {
    // Assing return storage
    BasisArray ret(mmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numgrid); // Potential
	u["density"  ].resize(numgrid); // Density
	u["rforce"   ].resize(numgrid); // Radial force
      }
    }

    // Radial grid spacing
    double dx = (logxmax - logxmin)/(numgrid-1);

    // Basis storage
    Eigen::MatrixXd tabpot, tabden, tabrfc;

    // Evaluate on the plane
    for (int i=0; i<numgrid; i++) {
      get_pot  (tabpot, pow(10.0, logxmin + dx*i));
      get_dens (tabden, pow(10.0, logxmin + dx*i));
      get_force(tabrfc, pow(10.0, logxmin + dx*i));
      for (int m=0; m<=mmax; m++) {
	for (int n=0; n<nmax; n++){
	  ret[m][n]["potential"](i) = tabpot(m, n);
	  ret[m][n]["density"  ](i) = tabden(m, n);
	  ret[m][n]["rforce"   ](i) = tabrfc(m, n);
	}
      }
    }
    
    return ret;
  }

  const std::set<std::string>
  Slab::valid_keys = {
    "nmaxx",
    "nmaxy",
    "nmaxz",
    "nminx",
    "nminy",
    "hslab",
    "zmax",
    "ngrid",
    "type",
    "knots",
    "verbose",
    "check",
    "method"
  };

  Slab::Slab(const YAML::Node& CONF) : BiorthBasis(CONF, "slab")
  {
    initialize();
  }

  Slab::Slab(const std::string& confstr) : BiorthBasis(confstr, "slab")
  {
    initialize();
  }

  void Slab::initialize()
  {
    // Basis identifier
    //
    BasisID = "Slab";

    nminx = 0;
    nminy = 0;

    nmaxx = 6;
    nmaxy = 6;
    nmaxz = 6;

    knots = 40;

    // Check orthogonality (false by default because of its long
    // runtime and very low utility)
    //
    bool check = false;

    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::Slab", "parameter", unmatched, __FILE__, __LINE__);
    
    // Default cachename, empty by default
    //
    std::string cachename;

    // Assign values from YAML
    //
    try {
      if (conf["nminx"])      nminx = conf["nminx"].as<int>();
      if (conf["nminy"])      nminy = conf["nminy"].as<int>();
      
      if (conf["nmaxx"])      nmaxx = conf["nmaxx"].as<int>();
      if (conf["nmaxy"])      nmaxy = conf["nmaxy"].as<int>();
      if (conf["nmaxz"])      nmaxz = conf["nmaxz"].as<int>();
      
      if (conf["hslab"])      hslab = conf["hslab"].as<double>();
      if (conf["zmax" ])      zmax  = conf["zmax" ].as<double>();
      if (conf["ngrid"])      ngrid = conf["ngrid"].as<int>();
      if (conf["type" ])      type  = conf["type" ].as<std::string>();

      if (conf["knots"])      knots = conf["knots"].as<int>();

      if (conf["check"])      check = conf["check"].as<bool>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Slab: error parsing YAML");
    }
    
    // Finally, make the basis
    //
    SLGridSlab::mpi  = 0;
    SLGridSlab::ZBEG = 0.0;
    SLGridSlab::ZEND = 0.1;
    SLGridSlab::H    = hslab;
  
    int nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

    ortho = std::make_shared<SLGridSlab>(nnmax, nmaxz, ngrid, zmax, type);

    // Orthogonality sanity check
    //
    if (check and myid==0) orthoTest();

    // Get max threads
    //
    int nthrds = omp_get_max_threads();

    imx = 2*nmaxx + 1;		// x wave numbers
    imy = 2*nmaxy + 1;		// y wave numbers
    imz = nmaxz;		// z basis count

    // Coefficient tensor
    //
    expcoef.resize(imx, imy, imz);
    expcoef.setZero();
      
    used = 0;

    // Set cartesian coordindates
    //
    coordinates = Coord::Cartesian;
  }
  
  void Slab::reset_coefs(void)
  {
    expcoef.setZero();
    totalMass = 0.0;
    used = 0;
  }
  
  
  void Slab::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    auto cf = dynamic_cast<CoefClasses::SlabStruct*>(coef.get());

    cf->nmaxx   = nmaxx;
    cf->nmaxy   = nmaxy;
    cf->nmaxz   = nmaxz;
    cf->time    = time;

    cf->allocate();

    *cf->coefs = expcoef;
  }

  void Slab::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (not dynamic_cast<CoefClasses::SlabStruct*>(coef.get()))
      throw std::runtime_error("Slab::set_coefs: you must pass a CoefClasses::SlabStruct");

    // Sanity check on dimensionality
    //
    {
      auto cc = dynamic_cast<CoefClasses::SlabStruct*>(coef.get());
      auto d  = cc->coefs->dimensions();
      if (d[0] != 2*nmaxx+1 or d[1] != 2*nmaxy+1 or d[2] != nmaxz) {
	std::ostringstream sout;
	sout << "Slab::set_coefs: the basis has (2*nmaxx+1, 2*nmaxy+1, nmaxz)=("
	     << 2*nmaxx+1 << ", " 
	     << 2*nmaxy+1 << ", " 
	     << nmaxz
	     << "). The coef structure has dimension=("
	     << d[0] << ", " << d[1] << ", " << d[2] << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    auto cf = dynamic_cast<CoefClasses::SlabStruct*>(coef.get());
    expcoef = *cf->coefs;

    // Cache the current coefficient structure
    //
    coefret = coef;

    coefctr = {0.0, 0.0, 0.0};
  }

  void Slab::accumulate(double x, double y, double z, double mass)
  {
    // Truncate to slab with sides in [0,1]
    if (x<0.0)
      x += std::floor(-x) + 1.0;
    else
      x -= std::floor( x);
    
    if (y<0.0)
      y += std::floor(-y) + 1.0;
    else
      y -= std::floor( y);
    
    used++;

    // Storage for basis evaluation
    Eigen::VectorXd zpot(nmaxz);

    // Loop indices
    int ix, iy;

    // Recursion multipliers
    std::complex<double> stepx = exp(-kfac*x), facx;
    std::complex<double> stepy = exp(-kfac*y), facy;
   
    // Initial values
    std::complex<double> startx = exp(static_cast<double>(nmaxx)*kfac*x);
    std::complex<double> starty = exp(static_cast<double>(nmaxy)*kfac*y);
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      // Wave number
      int ii = ix - nmaxx;
      int iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	// Wave number
	int jj = iy - nmaxy;
	int iiy = abs(jj);
	
	if (iix > nmaxx) {
	  std::cerr << "Out of bounds: iix=" << ii << std::endl;
	}
	if (iiy > nmaxy) {
	  std::cerr << "Out of bounds: iiy=" << jj << std::endl;
	}
	
	// Evaluate basis
	if (iix>=iiy)
	  ortho->get_pot(zpot, z, iix, iiy);
	else
	  ortho->get_pot(zpot, z, iiy, iix);

	for (int iz=0; iz<imz; iz++) {
	                       // +--- density in orthogonal series
                               // |    is 4.0*M_PI rho
                               // v
	  expcoef(ix, iy, iz) += -4.0*M_PI*mass*facx*facy*zpot[iz];
	}
      }
    }
  }
  
  void Slab::make_coefs()
  {
    if (use_mpi) {
      
      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      MPI_Allreduce(MPI_IN_PLACE, expcoef.data(), expcoef.size(), MPI_DOUBLE_COMPLEX,
		    MPI_SUM, MPI_COMM_WORLD);
    }
  }
  
  std::tuple<double, double, double, double, double>
  Slab::eval(double x, double y, double z)
  {
    // Loop indices
    //
    int ix, iy, iz;

    // Working values
    //
    std::complex<double> facx, facy, fac, facf, facd;

    // Return values
    //
    std::complex<double> accx(0.0), accy(0.0), accz(0.0), potl(0.0), dens(0.0);
    
    // Recursion multipliers
    //
    std::complex<double> stepx = exp(kfac*x);
    std::complex<double> stepy = exp(kfac*y);

    // Initial values (note sign change)
    //
    std::complex<double> startx = exp(-static_cast<double>(nmaxx)*kfac*x);
    std::complex<double> starty = exp(-static_cast<double>(nmaxy)*kfac*y);
    
    Eigen::VectorXd vpot(nmaxz), vfrc(nmaxz), vden(nmaxz);

    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      // Compute wavenumber; recall that the coefficients are stored
      // as follows: -nmax,-nmax+1,...,0,...,nmax-1,nmax
      //
      int ii = ix - nmaxx;
      int iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	int jj = iy - nmaxy;
	int iiy = abs(jj);
	
	if (iix > nmaxx) {
	  std::cerr << "Out of bounds: ii=" << ii << std::endl;
	}
	if (iiy > nmaxy) {
	  std::cerr << "Out of bounds: jj=" << jj << std::endl;
	}
	
	if (iix>=iiy) {
	  ortho->get_pot  (vpot, z, iix, iiy);
	  ortho->get_force(vfrc, z, iix, iiy);
	  ortho->get_dens (vden, z, iix, iiy);
	}
	else {
	  ortho->get_pot  (vpot, z, iiy, iix);
	  ortho->get_force(vfrc, z, iiy, iix);
	  ortho->get_dens (vden, z, iiy, iix);
	}

	
	for (int iz=0; iz<imz; iz++) {
	  
	  fac  = facx*facy*vpot[iz]*expcoef(ix, iy, iz);
	  facf = facx*facy*vfrc[iz]*expcoef(ix, iy, iz);
	  facd = facx*facy*vden[iz]*expcoef(ix, iy, iz);
	  
	  // Limit to minimum wave number
	  //
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  potl +=  fac;
	  dens +=  facd;
	  accx += -kfac*static_cast<double>(ii)*fac;
	  accy += -kfac*static_cast<double>(jj)*fac;
	  accz += -facf;
	  
	}
      }
    }

    return {potl.real(), dens.real(), accx.real(), accy.real(), accz.real()};
  }


  void Slab::computeAccel(double x, double y, double z,
			  Eigen::Ref<Eigen::Vector3d> acc)
  {
    // Loop indices
    //
    int ix, iy, iz;

    // Working values
    //
    std::complex<double> facx, facy, fac, facf;

    // Return values
    //
    std::complex<double> accx(0.0), accy(0.0), accz(0.0);
    
    // Recursion multipliers
    //
    std::complex<double> stepx = exp(kfac*x);
    std::complex<double> stepy = exp(kfac*y);

    // Initial values (note sign change)
    //
    std::complex<double> startx = exp(-static_cast<double>(nmaxx)*kfac*x);
    std::complex<double> starty = exp(-static_cast<double>(nmaxy)*kfac*y);
    
    Eigen::VectorXd vpot(nmaxz), vfrc(nmaxz);

    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      // Compute wavenumber; recall that the coefficients are stored
      // as follows: -nmax,-nmax+1,...,0,...,nmax-1,nmax
      //
      int ii = ix - nmaxx;
      int iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	int jj = iy - nmaxy;
	int iiy = abs(jj);
	
	if (iix > nmaxx) {
	  std::cerr << "Out of bounds: ii=" << ii << std::endl;
	}
	if (iiy > nmaxy) {
	  std::cerr << "Out of bounds: jj=" << jj << std::endl;
	}
	
	if (iix>=iiy) {
	  ortho->get_pot  (vpot, z, iix, iiy);
	  ortho->get_force(vfrc, z, iix, iiy);
	}
	else {
	  ortho->get_pot  (vpot, z, iiy, iix);
	  ortho->get_force(vfrc, z, iiy, iix);
	}

	
	for (int iz=0; iz<imz; iz++) {
	  
	  fac  = facx*facy*vpot[iz]*expcoef(ix, iy, iz);
	  facf = facx*facy*vfrc[iz]*expcoef(ix, iy, iz);
	  
	  // Limit to minimum wave number
	  //
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  accx += -kfac*static_cast<double>(ii)*fac;
	  accy += -kfac*static_cast<double>(jj)*fac;
	  accz += -facf;
	  
	}
      }
    }

    acc << accx.real(), accy.real(), accz.real();
  }


  std::vector<double> Slab::crt_eval(double x, double y, double z)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    auto [pot, den, frcx, frcy, frcz] = eval(x, y, z);

    return {0, den, den, 0, pot, pot, frcx, frcy, frcz};
  }

  std::vector<double> Slab::cyl_eval(double R, double z, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Cartesian from Cylindrical coordinates
    double x = R*cos(phi), y = R*sin(phi);

    auto [pot, den, frcx, frcy, frcz] = eval(x, y, z);

    double potR =  frcx*cos(phi) + frcy*sin(phi);
    double potp = -frcx*sin(phi) + frcy*cos(phi);
    double potz =  frcz;

    potR *= -1;
    potp *= -1;
    potz *= -1;

    return {0, den, den, 0, pot, pot, potR, potz, potp};
  }

  std::vector<double> Slab::sph_eval(double r, double costh, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Spherical from Cylindrical coordinates
    double sinth = sqrt(fabs(1.0 - costh*costh));
    double x = r*cos(phi)*sinth, y = r*sin(phi)*sinth, z = r*costh;

    auto [pot, den, frcx, frcy, frcz] = eval(x, y, z);

    double potr =  frcx*cos(phi)*sinth + frcy*sin(phi)*sinth + frcz*costh;
    double pott =  frcx*cos(phi)*costh + frcy*sin(phi)*costh - frcz*sinth;
    double potp = -frcx*sin(phi)       + frcy*cos(phi);

    potr *= -1;
    pott *= -1;
    potp *= -1;
    
    return {0, den, den, 0, pot, pot, potr, pott, potp};
  }


  Slab::BasisArray Slab::getBasis
  (double zmin, double zmax, int numgrid)
  {
    // Assign storage for returned basis array.  The Maximum
    // wavenumber for the SLGridSlab is the maximum of the X and Y
    // wavenumbers.
    int nnmax = std::max<int>(nmaxx, nmaxy);

    BasisArray ret (nnmax+1);	// X wavenumbers
    for (auto & v1 : ret) {
      v1.resize(nnmax+1);	// Y wavenumbers

      for (auto & v2 : v1) {
	v2.resize(nmaxz);	// Z basis

	for (auto & u : v2) {
	  u["potential"].resize(numgrid); // Potential
	  u["density"  ].resize(numgrid); // Density
	  u["zforce"   ].resize(numgrid); // Vertical force
	}
      }
    }

    // Vertical grid spacing
    double dz = (zmax - zmin)/(numgrid-1);

    // Basis evaluation storage
    Eigen::VectorXd vpot(nmaxz), vfrc(nmaxz), vden(nmaxz);

    // Construct the tensor
    for (int ix=0; ix<=nnmax; ix++) {
      
      for (int iy=0; iy<=nmaxx; iy++) {
      
	for (int i=0; i<numgrid; i++) {

	  double z = zmin + dz*i;

	  if (ix>=iy) {
	    ortho->get_pot  (vpot, z, ix, iy);
	    ortho->get_force(vfrc, z, ix, iy);
	    ortho->get_dens (vden, z, ix, iy);
	  } else {
	    ortho->get_pot  (vpot, z, iy, ix);
	    ortho->get_force(vfrc, z, iy, ix);
	    ortho->get_dens (vden, z, iy, ix);
	  }

	  for (int n=0; n<nmaxz; n++){
	    ret[ix][iy][n]["potential"](i) = vpot(n);
	    ret[ix][iy][n]["density"  ](i) = vden(n);
	    ret[ix][iy][n]["zforce"   ](i) = vfrc(n);
	  }
	}
      }
    }
    
    // Return the tensor
    return ret;
  }

  std::vector<Eigen::MatrixXd> Slab::orthoCheck(int knots)
  {
    return ortho->orthoCheck();
  }
  
  const std::set<std::string>
  Cube::valid_keys = {
    "nminx",
    "nminy",
    "nminz",
    "nmaxx",
    "nmaxy",
    "nmaxz",
    "knots",
    "verbose",
    "check",
    "method"
  };

  Cube::Cube(const YAML::Node& CONF) : BiorthBasis(CONF, "cube")
  {
    initialize();
  }

  Cube::Cube(const std::string& confstr) : BiorthBasis(confstr, "cube")
  {
    initialize();
  }

  void Cube::initialize()
  {
    // Basis identifier
    //
    BasisID = "Cube";

    nminx = std::numeric_limits<int>::max();
    nminy = std::numeric_limits<int>::max();
    nminz = std::numeric_limits<int>::max();

    nmaxx = 6;
    nmaxy = 6;
    nmaxz = 6;

    knots = 40;

    // Check orthogonality (false by default because of its long
    // runtime and very low utility)
    //
    bool check = false;

    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::Cube", "parameter", unmatched, __FILE__, __LINE__);
    
    // Default cachename, empty by default
    //
    std::string cachename;

    // Assign values from YAML
    //
    try {
      if (conf["nminx"])      nminx = conf["nminx"].as<int>();
      if (conf["nminy"])      nminy = conf["nminy"].as<int>();
      if (conf["nminz"])      nminz = conf["nminz"].as<int>();
      
      if (conf["nmaxx"])      nmaxx = conf["nmaxx"].as<int>();
      if (conf["nmaxy"])      nmaxy = conf["nmaxy"].as<int>();
      if (conf["nmaxz"])      nmaxz = conf["nmaxz"].as<int>();
      
      if (conf["knots"])      knots = conf["knots"].as<int>();

      if (conf["check"])      check = conf["check"].as<bool>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Cube: error parsing YAML");
    }
    
    // Finally, make the basis
    //
    ortho = std::make_shared<BiorthCube>(conf);
    
    // Orthogonality sanity check
    //
    if (check and myid==0) orthoTest();

    // Get max threads
    //
    int nthrds = omp_get_max_threads();

    expcoef.resize(2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1);
    expcoef.setZero();
      
    used = 0;

    // Set cartesian coordindates
    //
    coordinates = Coord::Cartesian;
  }
  
  void Cube::reset_coefs(void)
  {
    expcoef.setZero();
    totalMass = 0.0;
    used = 0;
  }
  
  
  void Cube::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    auto cf = dynamic_cast<CoefClasses::CubeStruct*>(coef.get());

    cf->nmaxx   = nmaxx;
    cf->nmaxy   = nmaxy;
    cf->nmaxz   = nmaxz;
    cf->time    = time;

    cf->allocate();

    *cf->coefs = expcoef;
  }

  void Cube::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (not dynamic_cast<CoefClasses::CubeStruct*>(coef.get()))
      throw std::runtime_error("Cube::set_coefs: you must pass a CoefClasses::CubeStruct");

    // Sanity check on dimensionality
    //
    {
      auto cc = dynamic_cast<CoefClasses::CubeStruct*>(coef.get());
      auto d  = cc->coefs->dimensions();
      if (d[0] != 2*nmaxx+1 or d[1] != 2*nmaxy+1 or d[2] != 2*nmaxz+1) {
	std::ostringstream sout;
	sout << "Cube::set_coefs: the basis has (2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1)=("
	     << 2*nmaxx+1 << ", " 
	     << 2*nmaxy+1 << ", " 
	     << 2*nmaxz+1
	     << "). The coef structure has dimension=("
	     << d[0] << ", " << d[1] << ", " << d[2] << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    auto cf = dynamic_cast<CoefClasses::CubeStruct*>(coef.get());
    expcoef = *cf->coefs;

    // Cache the cuurent coefficient structure
    //
    coefret = coef;

    coefctr = {0.0, 0.0, 0.0};
  }

  void Cube::accumulate(double x, double y, double z, double mass)
  {
    // Truncate to cube with sides in [0,1]
    if (x<0.0)
      x += std::floor(-x) + 1.0;
    else
      x -= std::floor( x);
    
    if (y<0.0)
      y += std::floor(-y) + 1.0;
    else
      y -= std::floor( y);
    
    if (z<0.0)
      z += std::floor(-z) + 1.0;
    else
      z -= std::floor( z);
    
    
    // Recursion multipliers
    Eigen::Vector3cd step
      {std::exp(-kfac*x), std::exp(-kfac*y), std::exp(-kfac*z)};
    
    // Initial values for recursion
    Eigen::Vector3cd init
      {std::exp(-kfac*(x*nmaxx)),
       std::exp(-kfac*(y*nmaxy)),
       std::exp(-kfac*(z*nmaxz))};
    
    Eigen::Vector3cd curr(init);
    for (int ix=0; ix<=2*nmaxx; ix++, curr(0)*=step(0)) {
      curr(1) = init(1);
      for (int iy=0; iy<=2*nmaxy; iy++, curr(1)*=step(1)) {
	curr(2) = init(2);
	for (int iz=0; iz<=2*nmaxz; iz++, curr(2)*=step(2)) {
	  
	  // Compute wavenumber; recall that the coefficients are
	  // stored as: -nmax,-nmax+1,...,0,...,nmax-1,nmax
	  //
	  int ii = ix-nmaxx;
	  int jj = iy-nmaxy;
	  int kk = iz-nmaxz;

	  // Normalization
	  double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));;

	  expcoef(ix, iy, iz) += - mass * curr(0)*curr(1)*curr(2) * norm;
	}
      }
    }
  }
  
  void Cube::make_coefs()
  {
    if (use_mpi) {
      
      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      MPI_Allreduce(MPI_IN_PLACE, expcoef.data(), expcoef.size(), MPI_DOUBLE_COMPLEX,
		    MPI_SUM, MPI_COMM_WORLD);
    }
  }
  
  std::vector<double> Cube::crt_eval(double x, double y, double z)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Position vector
    Eigen::Vector3d pos {x, y, z};

    // Get the basis fields
    double den1 = ortho->get_dens(expcoef, pos).real();
    double pot1 = ortho->get_pot (expcoef, pos).real();

    auto frc = ortho->get_force(expcoef, pos);
    
    double frcx = -frc(0).real();
    double frcy = -frc(1).real();
    double frcz = -frc(2).real();

    return {0, den1, den1, 0, pot1, pot1, frcx, frcy, frcz};
  }

  void Cube::computeAccel(double x, double y, double z,
			  Eigen::Ref<Eigen::Vector3d> acc)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Position vector
    Eigen::Vector3d pos {x, y, z};

    // Get the basis fields
    auto frc = ortho->get_force(expcoef, pos);
    
    acc << -frc(0).real(), -frc(1).real(), -frc(2).real();
  }

  std::vector<double> Cube::cyl_eval(double R, double z, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Cartesian from Cylindrical coordinates
    double x = R*cos(phi), y = R*sin(phi);

    // Position vector
    Eigen::Vector3d pos {x, y, z};

    // Get the basis fields
    double den1 = ortho->get_dens(expcoef, pos).real();
    double pot1 = ortho->get_pot (expcoef, pos).real();

    auto frc = ortho->get_force(expcoef, pos);
    
    double frcx = frc(0).real(), frcy = frc(1).real(), frcz = frc(2).real();

    double potR =  frcx*cos(phi) + frcy*sin(phi);
    double potp = -frcx*sin(phi) + frcy*cos(phi);
    double potz =  frcz;

    potR *= -1;
    potp *= -1;
    potz *= -1;

    return {0, den1, den1, 0, pot1, pot1, potR, potz, potp};
  }

  std::vector<double> Cube::sph_eval(double r, double costh, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Spherical from Cylindrical coordinates
    double sinth = sqrt(fabs(1.0 - costh*costh));
    double x = r*cos(phi)*sinth, y = r*sin(phi)*sinth, z = r*costh;

    // Position vector
    Eigen::Vector3d pos {x, y, z};

    // Get the basis fields
    double den1 = ortho->get_dens(expcoef, pos).real();
    double pot1 = ortho->get_pot (expcoef, pos).real();

    auto frc = ortho->get_force(expcoef, pos);
    
    double frcx = frc(0).real();
    double frcy = frc(1).real();
    double frcz = frc(2).real();

    double potr =  frcx*cos(phi)*sinth + frcy*sin(phi)*sinth + frcz*costh;
    double pott =  frcx*cos(phi)*costh + frcy*sin(phi)*costh - frcz*sinth;
    double potp = -frcx*sin(phi)       + frcy*cos(phi);

    potr *= -1;
    pott *= -1;
    potp *= -1;
    
    return {0, den1, den1, 0, pot1, pot1, potr, pott, potp};
  }

  std::vector<Eigen::MatrixXd> Cube::orthoCheck(int knots)
  {
    std::vector<Eigen::MatrixXd> ret;
    ret.push_back(ortho->orthoCheck().array().abs());
    return ret;
  }
  
  // Generate coeffients from a particle reader
  CoefClasses::CoefStrPtr BiorthBasis::createFromReader
  (PR::PRptr reader, Eigen::Vector3d ctr, RowMatrix3d rot)
  {
    CoefClasses::CoefStrPtr coef;

    if (name.compare("sphereSL") == 0)
      coef = std::make_shared<CoefClasses::SphStruct>();
    else if (name.compare("cylinder") == 0)
      coef = std::make_shared<CoefClasses::CylStruct>();
    else if (name.compare("flatdisk") == 0)
      coef = std::make_shared<CoefClasses::CylStruct>();
    else if (name.compare("cube") == 0)
      coef = std::make_shared<CoefClasses::CubeStruct>();
    else {
      std::ostringstream sout;
      sout << "Basis::createCoefficients: basis <" << name << "> not recognized"
	   << std::endl;
      throw std::runtime_error(sout.str());
    }
      
    // Add the expansion center metadata and register for this instance
    //
    coefctr   = ctr;
    coef->ctr = ctr;

    // Add the rotation matrix metadata and register for this instance
    //
    coefrot   = rot;
    coef->rot = rot;

    std::vector<double> p1(3), v1(3);
    
    // Map the vector rather than copy
    //
    Eigen::Map<Eigen::Vector3d> pp(p1.data(), 3), vv(v1.data(), 3);
    vv.setZero();

    reset_coefs();
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {

      bool use = false;
      
      // Translate and rotate the position vector
      //
      for (int k=0; k<3; k++) p1[k] = p->pos[k];
      pp = coefrot * (pp - coefctr);

      if (ftor) {
	// Rotate the velocity vector
	//
	for (int k=0; k<3; k++) v1[k] = p->vel[k];
	vv = coefrot * vv;
	
	use = ftor(p->mass, p1, v1, p->indx);
      } else {
	use = true;
      }

      if (use) {
	accumulate(pp(0), pp(1), pp(2), p->mass);
      }
    }
    make_coefs();
    load_coefs(coef, reader->CurrentTime());
    return coef;
  }

  // Generate coefficients from a phase-space table
  void BiorthBasis::initFromArray(Eigen::Vector3d ctr, RowMatrix3d rot)
  {
    if (name.compare("sphereSL") == 0)
      coefret = std::make_shared<CoefClasses::SphStruct>();
    else if (name.compare("cylinder") == 0)
      coefret = std::make_shared<CoefClasses::CylStruct>();
    else if (name.compare("flatdisk") == 0)
      coefret = std::make_shared<CoefClasses::CylStruct>();
    else {
      std::ostringstream sout;
      sout << "Basis::createCoefficients: basis <" << name << "> not recognized"
	   << std::endl;
      throw std::runtime_error(sout.str());
    }
      
    // Add the expansion center metadata and register
    //
    coefctr = ctr;
    coefret->ctr = ctr;

    // Add the rotation metadata and register
    coefrot = rot;
    coefret->rot = rot;

    // Clean up for accumulation
    //
    reset_coefs();
    coefindx = 0;
  }

  // Accumulate coefficient contributions from arrays
  void BiorthBasis::addFromArray(Eigen::VectorXd& m, RowMatrixXd& p,
				 bool RoundRobin, bool PosVelRows)
  {
    // Sanity check: is coefficient instance created?  This is not
    // foolproof.  It is really up the user to make sure that a call
    // to initFromArray() comes first.
    //
    if (not coefret) {
      std::string msg =
	"Basis::addFromArray: you must initialize coefficient accumulation "
	"with a call to Basis::initFromArray()";
      throw std::runtime_error(msg);
    }

    // Assume position arrays in rows by default
    //
    int rows = p.rows();
    int cols = p.cols(); 

    bool ambiguous = false;
    bool haveVel   = false;

    if (cols==3 or cols==6) {
      if (rows != 3 and rows != 6) PosVelRows = false;
      else ambiguous = true;
    }

    if (rows==3 or rows==6) {
      if (cols != 3 and cols != 6) PosVelRows = true;
      else ambiguous = true;
    }

    if (ambiguous and myid==0) {
      std::cout << "---- BiorthBasis::addFromArray: dimension deduction "
		<< "is ambiguous.  Assuming that ";
      if (PosVelRows) std::cout << "positions are in rows" << std::endl;
      else std::cout << "positions are in columns" << std::endl;
      std::cout << "---- BiorthBasis::addFromArray: reset 'posvelrows' flag "
		<< "if this assumption is wrong." << std::endl;
    }

    // Map the vector rather than copy
    //
    std::vector<double> p1(3), v1(3);
    Eigen::Map<Eigen::Vector3d> pp(p1.data(), 3), vv(v1.data(), 3);
    vv.setZero();

    if (PosVelRows) {
      if (p.rows()<3) {
	std::ostringstream msg;
	msg << "Basis::addFromArray: you must pass a position array with at "
	  "least three rows for x, y, z.  Yours has " << p.rows() << ".";
	throw std::runtime_error(msg.str());
      }

      if (p.rows() == 6) haveVel = true;

      for (int n=0; n<p.cols(); n++) {
	
	if (n % numprocs==myid or not RoundRobin) {

	  for (int k=0; k<3; k++) {
	    pp(k) = p(k, n);
	    if (haveVel) vv(k) = p(k+3, n);
	  }
	  
	  pp = coefrot * (pp - coefctr);

	  bool use = true;

	  if (ftor) {
	    if (haveVel) vv = coefrot * vv;
	    use = ftor(m(n), p1, v1, coefindx);
	  } else {
	    use = true;
	  }
	  coefindx++;
	  
	  if (use) {
	    accumulate(pp(0), pp(1), pp(2), m(n));
	  }
	}
      }
      
    } else {

      if (p.cols()<3) {
	std::ostringstream msg;
	msg << "Basis::addFromArray: you must pass a position array with at "
	  "least three columns for x, y, z.  Yours has " << p.cols() << ".";
	throw std::runtime_error(msg.str());
      }

      if (p.cols() == 6) haveVel = true;


      for (int n=0; n<p.rows(); n++) {

	if (n % numprocs==myid or not RoundRobin) {

	  for (int k=0; k<3; k++) {
	    pp(k) = p(n, k);
	    if (haveVel) vv(k) = p(n, k+3);
	  }
	  
	  pp = coefrot * (pp - coefctr);

	  bool use = true;
	  if (ftor) {
	    if (haveVel) vv = coefrot * vv;
	    use = ftor(m(n), p1, v1, coefindx);
	  } else {
	    use = true;
	  }
	  coefindx++;
	  
	  if (use) {
	    accumulate(pp(0), pp(1), pp(2), m(n));
	  }
	}
      }
    }
  }

  // Generate coefficients from the accumulated array values
  CoefClasses::CoefStrPtr BiorthBasis::makeFromArray(double time)
  {
    make_coefs();
    load_coefs(coefret, time);
    return coefret;
  }

  // Generate coefficients from a phase-space table
  //
  CoefClasses::CoefStrPtr BiorthBasis::createFromArray
  (Eigen::VectorXd& m, RowMatrixXd& p, double time, Eigen::Vector3d ctr,
   RowMatrix3d rot, bool RoundRobin, bool PosVelRows)
  {
    initFromArray(ctr, rot);
    addFromArray(m, p, RoundRobin, PosVelRows);
    return makeFromArray(time);
  }

  // This evaluation step is performed by all derived classes
  Eigen::MatrixXd& AccelFunc::evalaccel
  (Eigen::MatrixXd& ps, Eigen::MatrixXd& accel, BasisCoef mod)
  {
    // Get Model
    //
    auto basis = std::get<0>(mod);

    // Get expansion center
    //
    auto ctr = basis->getCenter();
    if (basis->usingNonInertial()) ctr = {0, 0, 0};

    // Get rotation matrix
    //
    auto rot = basis->getRotation();

    // Get fields
    //
    int rows = accel.rows();
    for (int n=0; n<rows; n++) {
      Eigen::Vector3d pp;
      for (int k=0; k<3; k++) pp(k) = ps(n, k) - ctr(k);
      pp = rot * pp;

      auto v = basis->getFields(pp(0), pp(1), pp(2));

      // First 6 fields are density and potential, followed by acceleration
      for (int k=0; k<3; k++) accel(n, k) += v[6+k] - basis->pseudo(k);
    }

    // true for deep debugging
    //  |
    //  v
    if (false and basis->usingNonInertial()) {

      auto coefs = basis->getCoefficients();
      auto time  = coefs->time;
      auto ctr   = coefs->ctr;

      std::ofstream tmp;
      if (time <= 0.0) tmp.open("pseudo.dat");
      else             tmp.open("pseudo.dat", ios::app);

      if (tmp)
	tmp << std::setw(16) << std::setprecision(5) << time
	    << std::setw(16) << std::setprecision(5) << ctr[0]
	    << std::setw(16) << std::setprecision(5) << ctr[1]
	    << std::setw(16) << std::setprecision(5) << ctr[2]
	    << std::setw(16) << std::setprecision(5) << basis->pseudo(0)
	    << std::setw(16) << std::setprecision(5) << basis->pseudo(1)
	    << std::setw(16) << std::setprecision(5) << basis->pseudo(2)
	    << std::endl;
    }

    return accel;
  }

  // This is an example of a evalcoefs() derived class.  It is the
  // responsibility of the derived-class implementer to provide a sane
  // set of coefficients using Basis::set_coefs for each
  // component. Although not needed here, the best way of identifying
  // the component might be to use the getName() member of coefs,
  // e.g. 'std::string name = std::get<mod>(1)->getName();'
  void
  AllTimeAccel::evalcoefs(double t, BasisCoef mod)
  {
    auto basis = std::get<0>(mod);
    auto coefs = std::get<1>(mod);

    // Interpolate coefficients
    //
    auto times = coefs->Times();

    if (t<times.front() or t>times.back()) {
      std::ostringstream sout;
      sout << "Basis::OneAccel: time t=" << t << " is out of bounds: ["
	   << times.front() << ", " << times.back() << "]";
      throw std::runtime_error(sout.str());
    }
    
    auto it2 = std::lower_bound(times.begin(), times.end(), t);
    auto it1 = it2;

    if (it2 == times.end()) throw std::runtime_error("Basis::AllTimeAccel::evalcoefs: time t=" + std::to_string(t) + " out of bounds");
    else if (it2 == times.begin()) it2++;
    else it1--;

    double a = (*it2 - t)/(*it2 - *it1);
    double b = (t - *it1)/(*it2 - *it1);

    auto coefsA = coefs->getCoefStruct(*it1);
    auto coefsB = coefs->getCoefStruct(*it2);

    // Duplicate a coefficient instance
    //
    auto newcoef = coefsA->deepcopy();

    // Now interpolate the matrix
    //
    newcoef->time = t;

    auto & cN = newcoef->store;
    auto & cA = coefsA->store;
    auto & cB = coefsB->store;

    for (int i=0; i<newcoef->store.size(); i++)
      cN(i) = a * cA(i) + b * cB(i);

    // Interpolate center
    //
    newcoef->ctr = a * coefsA->ctr + b * coefsB->ctr;

    // Interpolate rotation matrix followed by unitarization
    //
    RowMatrix3d newrot = a * coefsA->rot + b * coefsB->rot;

    // Closest unitary matrix in the Frobenius norm sense
    //
    Eigen::BDCSVD<RowMatrix3d> svd
      (newrot, Eigen::ComputeFullU | Eigen::ComputeFullV);

    newcoef->rot = svd.matrixU() * svd.matrixV().adjoint();

    // Install coefficients
    //
    basis->set_coefs(newcoef);

    // Set non-inertial force
    basis->setNonInertialAccel(t);

  }

  SingleTimeAccel::SingleTimeAccel(double t, std::vector<BasisCoef> mod)
  {
    for (auto model : mod) {

      auto basis = std::get<0>(model);
      auto coefs = std::get<1>(model);

      // Interpolate coefficients
      //
      auto times = coefs->Times();

      if (t<times.front() or t>times.back()) {
	std::ostringstream sout;
	sout << "Basis::OneAccel: time t=" << t << " is out of bounds: ["
	     << times.front() << ", " << times.back() << "]";
	throw std::runtime_error(sout.str());
      }
      
      auto it2 = std::lower_bound(times.begin(), times.end(), t);
      auto it1 = it2;

      if (it2 == times.end())
	throw std::runtime_error("Basis::SingleTimeAccel::evalcoefs: time t=" + std::to_string(t) + " out of bounds");
      else if (it2 == times.begin()) it2++;
      else it1--;
      
      double a = (*it2 - t)/(*it2 - *it1);
      double b = (t - *it1)/(*it2 - *it1);
      
      auto coefsA = coefs->getCoefStruct(*it1);
      auto coefsB = coefs->getCoefStruct(*it2);
      
      // Duplicate a coefficient instance
      //
      auto newcoef = coefsA->deepcopy();
      
      // Now interpolate the matrix
      //
      newcoef->time = t;

      auto & cN = newcoef->store;
      auto & cA = coefsA ->store;
      auto & cB = coefsB ->store;

      for (int i=0; i<cN.size(); i++)
	cN(i) = a * cA(i) + b * cB(i);

      // Interpolate center
      //
      if (coefsA->ctr.size() and coefsB->ctr.size()) {
	newcoef->ctr.resize(3);
	for (int k=0; k<3; k++)
	  newcoef->ctr[k] = a * coefsA->ctr[k] + b * coefsB->ctr[k];
      }

      // Install coefficients
      //
      basis->set_coefs(newcoef);
    }
    // END: component model loop
  }
  
  //! Take one leap frog step; this can/should be generalized to a
  //! one-step class in the long run
  std::tuple<double, Eigen::MatrixXd>
  OneStep(double t, double h,
	  Eigen::MatrixXd& ps, Eigen::MatrixXd& accel,
	  std::vector<BasisCoef> bfe, AccelFunctor F)
  {
    int rows = ps.rows();

    // Leap frog (set to false for RK4 test)
    //
    if (true) {

      // Drift 1/2
      for (int n=0; n<rows; n++) {
	for (int k=0; k<3; k++) ps(n, k) += ps(n, 3+k)*0.5*h;
      }

      // Kick
      accel.setZero();
      for (auto mod : bfe) F(t, ps, accel, mod);

      for (int n=0; n<rows; n++) {
	for (int k=0; k<3; k++) ps(n, 3+k) += accel(n, k)*h;
      }
      
      // Drift 1/2
      for (int n=0; n<rows; n++) {
	for (int k=0; k<3; k++) ps(n, k) += ps(n, 3+k)*0.5*h;
      }
    }
    // RK4
    else {
      // Make and clear variables
      std::vector<Eigen::MatrixXd> kf(4);
      for (int i=0; i<4; i++) kf[i].resize(rows, 6);

      // Step 1
      //
      accel.setZero();
      for (auto mod : bfe) F(t, ps, accel, mod); 
      for (int n=0; n<rows; n++) {
	for (int k=0; k<3; k++) {
	  kf[0](n, 0+k) = ps(n, 3+k);
	  kf[0](n, 3+k) = accel(n, k);
	}
      }

      // Step 2
      //
      Eigen::MatrixXd ps1 = ps + kf[0]*0.5*h; // state vector update

      accel.setZero();
      for (auto mod : bfe) F(t+0.5*h, ps1, accel, mod); 
      for (int n=0; n<rows; n++) {
	for (int k=0; k<3; k++) {
	  kf[1](n, 0+k) = ps1(n, 3+k);
	  kf[1](n, 3+k) = accel(n, k);
	}
      }

      // Step 3
      //
      ps1 = ps + kf[1]*0.5*h;	// state vector update

      accel.setZero();
      for (auto mod : bfe) F(t+0.5*h, ps1, accel, mod); 
      for (int n=0; n<rows; n++) {
	for (int k=0; k<3; k++) {
	  kf[2](n, 0+k) = ps1(n, 3+k);
	  kf[2](n, 3+k) = accel(n, k);
	}
      }

      // Step 4
      ps1 = ps + kf[2]*h;	// state vector update

      accel.setZero();
      for (auto mod : bfe) F(t+h, ps1, accel, mod); 
      for (int n=0; n<rows; n++) {
	for (int k=0; k<3; k++) {
	  kf[3](n, 0+k) = ps1(n, 3+k);
	  kf[3](n, 3+k) = accel(n, k);
	}
      }

      // The final lhs
      //
      Eigen::MatrixXd acc = (kf[0] + 2.0*kf[1] + 2.0*kf[2] + kf[3])/6.0;

      for (int n=0; n<rows; n++) { // Copy back acceleration for this step
	for (int k=0; k<3; k++) accel(n, k) = acc(n, 3+k);
      }

      ps += acc*h;
    }

    return std::tuple<double, Eigen::MatrixXd>(t+h, ps);
  }


  std::tuple<Eigen::VectorXd, Eigen::Tensor<float, 3>>
  IntegrateOrbits
  (double tinit, double tfinal, double h,
   Eigen::MatrixXd ps, std::vector<BasisCoef> bfe, AccelFunctor F,
   int nout)
  {
    int rows = ps.rows();
    int cols = ps.cols();

    // ps should be a (n, 6) table of phase-space initial conditions
    //
    if (cols != 6) {
      std::ostringstream sout;
      sout << "IntegrateOrbits: phase space array should be n x 6 where n is "
	   << "the number of particles.  You specified " << cols << " columns";
      throw std::runtime_error(sout.str());
    }

    // Allocate the acceleration array
    //
    Eigen::MatrixXd accel(rows, 3);

    // Sanity check
    //
    if (tfinal == tinit) {
      throw std::runtime_error
	("BasisClasses::IntegrateOrbits: tinit cannot be equal to tfinal");
    }

    if (h < 0.0 and tfinal > tinit) {
      throw std::runtime_error
	("BasisClasses::IntegrateOrbits: tfinal must be smaller than tinit "
	 "when step size is negative");
    }

    if (h > 0.0 and tfinal < tinit) {
      throw std::runtime_error
	("BasisClasses::IntegrateOrbits: tfinal must be larger than "
	 "tinit when step size is positive");
    }

    if ( (tfinal - tinit)/h >
	 static_cast<double>(std::numeric_limits<int>::max()) )
      {
	std::cout << "BasisClasses::IntegrateOrbits: step size is too small or "
		  << "time interval is too large." << std::endl;
	// Return empty data
	//
	return {Eigen::VectorXd(), Eigen::Tensor<float, 3>()};
      }
    
    // Number of steps
    //
    int numT = std::ceil( (tfinal - tinit)/h + 0.5);

    // Want both end points in the output at minimum
    //
    numT = std::max(2, numT);

    // Number of output steps
    //
    int stride = 1;		// Default stride
    if (nout>0) {		// User has specified output count...
      nout = std::max(2, nout);
      stride = std::ceil(static_cast<double>(numT)/static_cast<double>(nout));
      numT = (nout-1) * stride + 1;
    } else {			// Otherwise, use the default output number
      nout = numT;		// with the default stride
    }

    // Compute the interval-matching step
    //
    h = (tfinal - tinit)/(numT-1);

    // DEBUG
    if (false) 
      std::cout << "BasisClasses::IntegrateOrbits: choosing nout=" << nout
		<< " numT=" << numT << " h=" << h << " stride=" << stride
		<< std::endl;

    // Return data
    //
    Eigen::Tensor<float, 3> ret;

    try {
      ret.resize(rows, 6, nout);
    }
    catch (const std::bad_alloc& e) {
      std::cout << "BasisClasses::IntegrateOrbits: memory allocation failed: "
		<< e.what() << std::endl
		<< "Your requested number of orbits and time steps requires "
		<< std::floor(4.0*rows*6*nout/1e9)+1 << " GB free memory"
		<< std::endl;

      // Return empty data
      //
      return {Eigen::VectorXd(), Eigen::Tensor<float, 3>()};
    }

    // Time array
    //
    Eigen::VectorXd times(nout);
    
    // Assign the initial point
    //
    times(0) = tinit;
    for (int n=0; n<rows; n++)
      for (int k=0; k<6; k++) ret(n, k, 0) = ps(n, k);

    // Sign of h
    int sgn = (0 < h) - (h < 0);

    // Set the counters
    double tnow = tinit;
    int s = 0, cnt = 1;

    // Do the integration using stride for output
    while (s++ < numT) {
      if ( (tfinal - tnow)*sgn < h*sgn) h = tfinal - tnow;
      std::tie(tnow, ps) = OneStep(tnow, h, ps, accel, bfe, F);
      if (cnt < nout and s % stride == 0) {
	times(cnt) = tnow;
	for (int n=0; n<rows; n++)
	  for (int k=0; k<6; k++) ret(n, k, cnt) = ps(n, k);
	cnt += 1;
      }
    }

    // Corrects round off at end point
    //
    times(nout-1) = tnow;
    for (int n=0; n<rows; n++)
      for (int k=0; k<6; k++) ret(n, k, nout-1) = ps(n, k);

    return {times, ret};
  }

  void Spherical::writeCovarH5Params(HighFive::File& file)
  {
    file.createAttribute<int>("lmax", HighFive::DataSpace::From(lmax)).write(lmax);
    file.createAttribute<int>("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
    file.createAttribute<double>("scale", HighFive::DataSpace::From(scale)).write(scale);
    file.createAttribute<double>("rmin", HighFive::DataSpace::From(rmin)).write(rmin);
    file.createAttribute<double>("rmax", HighFive::DataSpace::From(rmax)).write(rmax);
  }
  
  void Cylindrical::writeCovarH5Params(HighFive::File& file)
  {
    file.createAttribute<int>("mmax", HighFive::DataSpace::From(mmax)).write(mmax);
    file.createAttribute<int>("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
    file.createAttribute<double>("scale", HighFive::DataSpace::From(rscl)).write(rscl);
    file.createAttribute<double>("rmin", HighFive::DataSpace::From(rmin)).write(rmin);
    file.createAttribute<double>("rmax", HighFive::DataSpace::From(rmax)).write(rmax);
    file.createAttribute<double>("acyl", HighFive::DataSpace::From(acyl)).write(acyl);
    file.createAttribute<double>("hcyl", HighFive::DataSpace::From(hcyl)).write(hcyl);
  }
  
  unsigned BiorthBasis::writeCovarH5(HighFive::Group& snaps, unsigned count, double time)
  {
    std::ostringstream stim;
    stim << std::setw(8) << std::setfill('0') << std::right << count++;
      
    // Make a new group for this time
    //
    HighFive::Group stanza = snaps.createGroup(stim.str());
      
    // Add a time attribute
    //
    stanza.createAttribute<double>("Time", HighFive::DataSpace::From(time)).write(time);
      
    // Add the sample statistics
    //
    HighFive::DataSet s1data = stanza.createDataSet("sampleCounts", sampleCounts);
    HighFive::DataSet s2data = stanza.createDataSet("sampleMasses", sampleMasses);

    auto covar = getCoefCovariance();

    // Add the samples
    //
    for (size_t T=0; T<sampleCounts.size(); T++) {

      // Group name
      std::ostringstream sT;
      sT << std::setw(8) << std::setfill('0') << std::right << T;
      
      // Make a new group
      HighFive::Group sample = stanza.createGroup(sT.str());

      // Pack the coefficient data
      size_t lmax = covar[T].size();
      size_t nmax = std::get<0>(covar[T][0]).rows();

      Eigen::VectorXd data(lmax*nmax);
      for (size_t l=0, c=0; l<lmax; l++) {
	for (size_t n=0; n<nmax; n++, c++) {
	  data(c) = std::get<0>(covar[T][l])(n);
	}
      }

      HighFive::DataSet coefdata = sample.createDataSet("coefficients", data);

      // Pack the covariance data in an upper triangular format
      //
      size_t diagonalSize = nmax*(nmax + 1)/2;
      data.resize(lmax*diagonalSize);
      for (size_t l=0, c=0; l<lmax; l++) {
	for (size_t n1=0; n1<nmax; n1++) {
	  for (size_t n2=n1; n2<nmax; n2++, c++) {
	    data(c) = std::get<1>(covar[T][l])(n1, n2);
	  }
	}
      }

      HighFive::DataSet covdata = sample.createDataSet("covariance", data);
    }
    // END: sample loop

    return count;
  }
  
  void BiorthBasis::writeCoefCovariance(const std::string& compname, const std::string& runtag, double time)
  {
    // Check check that variance computation is on
    //
    if (not pcavar) {
      std::cout << "BiorthBasis::writeCoefCovariance: covariance computation is disabled.  "
		<< "Set 'pcavar: true' to enable." << std::endl;
      return;
    }

    // Only root process writes
    //
    if (myid) return;

    // Check that there is something to write
    //
    int totalCount = 0;
    for (auto c : sampleCounts) totalCount += c;
    if (totalCount==0) {
      std::cout << "BiorthBasis::writeCoefCovariance: no data" << std::endl;
      return;
    }

    // The H5 filename
    std::string fname = "coefcovar." + compname + "." + runtag + ".h5";

    // Check if file exists?
    //
    try {
      // Open the HDF5 file in read-only mode
      HighFive::File file(fname, HighFive::File::ReadOnly);

      // Check for version string
      std::string path = "CovarianceFileVersion"; 

      // Check for valid HDF file by attribute
      if (file.hasAttribute(path)) {
	extendCoefCovariance(fname, time);
	return;
      }

    } catch (const HighFive::Exception& err) {
        // Handle HighFive specific errors (e.g., file not found)
        std::cerr << "HighFive Error: " << err.what() << std::endl;
    } catch (const std::exception& err) {
        // Handle other general exceptions
        std::cerr << "Error: " << err.what() << std::endl;
    }

    // Write coefficients
    //
    try {
      // Create a new hdf5 file
      //
      HighFive::File file(fname,
			  HighFive::File::ReadWrite |
			  HighFive::File::Create);
      
      // Write the Version string
      //
      file.createAttribute<std::string>("CovarianceFileVersion", HighFive::DataSpace::From(CovarianceFileVersion)).write(CovarianceFileVersion);

      // We write the basis identifier string
      //
      file.createAttribute<std::string>("BasisID", HighFive::DataSpace::From(BasisID)).write(BasisID);
      
      // Write the specific parameters
      //
      writeCovarH5Params(file);
      
      // Group count variable
      //
      unsigned count = 0;
      HighFive::DataSet dataset = file.createDataSet("count", count);
      
      // Create a new group for coefficient snapshots
      //
      HighFive::Group group = file.createGroup("snapshots");
      
      // Write the coefficients
      //
      count = writeCovarH5(group, count, time);
      
      // Update the count
      //
      dataset.write(count);
      
    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
  }
  
  void BiorthBasis::extendCoefCovariance(const std::string& fname, double time)
  {
    try {
      // Open an hdf5 file
      //
      HighFive::File file(fname, HighFive::File::ReadWrite);
      
      // Get the dataset
      HighFive::DataSet dataset = file.getDataSet("count");
      
      unsigned count;
      dataset.read(count);
      
      HighFive::Group group = file.getGroup("snapshots");
      
      // Write the coefficients
      //
      count = writeCovarH5(group, count, time);
      
      // Update the count
      //
      dataset.write(count);
      
    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
  }

  // Read covariance data
  CovarianceReader::CovarianceReader(const std::string& filename, int stride)
  {
    try {
      // Open an existing hdf5 file for reading
      //
      HighFive::File file(filename, HighFive::File::ReadOnly);
      
      // Write the Version string
      //
      std::string version;
      file.getAttribute("CovarianceFileVersion").read(version);
      if (version != std::string("1.0")) {
	throw std::runtime_error("CovarianceReader: unsupported file version, "
				 + version);
      }

      // Read the basis identifier string
      //
      file.getAttribute("BasisID").read(basisID);
      
      int lmax, nmax, ltot;

      // Current implemented spherical types
      const std::set<std::string> sphereType = {"Spherical", "SphereSL", "bessel"};

      // Currently implemented cylindrical types
      const std::set<std::string> cylinderType = {"Cylindrical"};

      if (sphereType.find(basisID) != sphereType.end()) {
	file.getAttribute("lmax").read(lmax);
	file.getAttribute("nmax").read(nmax);
	ltot = (lmax+1)*(lmax+2)/2;
      } else if (cylinderType.find(basisID) != cylinderType.end()) {
	file.getAttribute("mmax").read(lmax);
	file.getAttribute("nmax").read(nmax);
	ltot = lmax + 1;
      } else {
	throw std::runtime_error("CovarianceReader: unknown or unimplemented covariance for basis type, " + basisID);
      }

      // Group count variable
      //
      unsigned count = 0;
      file.getDataSet("count").read(count);

      // Open the snapshot group
      //
      auto snaps = file.getGroup("snapshots");
      
      for (unsigned n=0; n<count; n+=stride) {

	std::ostringstream sout;
	sout << std::setw(8) << std::setfill('0') << std::right << n;
      
	auto stanza = snaps.getGroup(sout.str());
      
	double Time;
	stanza.getAttribute("Time").read(Time);

	int itime = static_cast<int>(Time * fixedPointPrecision + 0.5);
	timeMap[itime] = times.size();
	times.push_back(Time);

	// Get sample properties
	//
	sampleCounts.push_back(Eigen::VectorXi());
	stanza.getDataSet("sampleCounts").read(sampleCounts.back());
	
	sampleMasses.push_back(Eigen::VectorXd());
	stanza.getDataSet("sampleMasses").read(sampleMasses.back());

	int nT = sampleCounts.back().size();

	// Allocate vector for current time
	covarData.push_back(std::vector<std::vector<CoefCovarType>>(nT));

	for (int T=0; T<nT; T++) {
	  // Group name
	  std::ostringstream sT;
	  sT << std::setw(8) << std::setfill('0') << std::right << T;
      
	  // Get the group
	  HighFive::Group sample = stanza.getGroup(sT.str());

	  // Storage
	  Eigen::VectorXd data;

	  // Repack the data
	  std::vector<CoefCovarType> elem(ltot);
	  for (auto & e : elem) {
	    std::get<0>(e).resize(nmax);
	    std::get<1>(e).resize(nmax, nmax);
	  }

	  // Get the flattened coefficient array
	  data = sample.getDataSet("coefficients").read<Eigen::VectorXd>();
	  
	  // Pack the coefficient data
	  for (size_t l=0, c=0; l<ltot; l++) {
	    for (size_t n=0; n<nmax; n++, c++) {
	      std::get<0>(elem[l])(n) = data(c);
	    }
	  }

	  // Get the flattened covariance array
	  data = sample.getDataSet("covariance").read<Eigen::VectorXd>();
	  
	  // Pack the coefficient data
	  for (size_t l=0, c=0; l<ltot; l++) {
	    for (size_t n1=0; n1<nmax; n1++) {
	      for (size_t n2=n1; n2<nmax; n2++, c++) {
		std::get<1>(elem[l])(n1, n2) = data(c);
		if (n1 != n2)
		  std::get<1>(elem[l])(n2, n1) = std::get<1>(elem[l])(n1, n2);
	      }
	    }
	  }

	  // Add the data
	  covarData.back()[T] = std::move(elem);
	}
	// END: sample loop

      }
      // END: snapshot loop
    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }
  }

}
// END namespace BasisClasses
