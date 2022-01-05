 /*
  Check orthgonality for cylindrical basis
*/
                                // C++/STL headers
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <cmath>

#include <fenv.h>

#include <yaml-cpp/yaml.h>	// YAML support

#include <config.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

// EXP classes
//
#include <numerical.H>
#include <gaussQ.H>
#include <isothermal.H>
#include <hernquist.H>
#include <model3d.H>
#include <biorth.H>
#include <SphericalSL.H>
#include <interp.H>
#include <EmpCylSL.H>
#include <DiskModels.H>
#include <DiskEval.H>
#include <norminv.H>
#include <cxxopts.H>
#include <EXPini.H>

#define M_SQRT1_3 (0.5773502691896257645091487)

                                // For debugging
#ifdef DEBUG
#include <fenv.h>
#include <fpetrap.h>

//===========================================
// Handlers defined in exputil/stack.cc
//===========================================

extern void mpi_print_trace(const string& routine, const string& msg,
			    const char *file, int line);

extern void mpi_gdb_print_trace(int sig);

extern void mpi_gdb_wait_trace(int sig);

//===========================================
// A signal handler to trap invalid FP only
//===========================================

void set_fpu_invalid_handler(void)
{
  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);
}

//===========================================
// A signal handler to produce a traceback
//===========================================

void set_fpu_trace_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);
}

//===========================================
// A signal handler to produce stop and wait
//===========================================

void set_fpu_gdb_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_wait_trace);
}

#endif

                                // Local headers
#include <localmpi.H>


// Hydrogen fraction
//
const double f_H = 0.76;


// Global variables
//
enum DiskType { constant, gaussian, mn, exponential };

std::map<std::string, DiskType> dtlookup =
  { {"constant",    DiskType::constant},
    {"gaussian",    DiskType::gaussian},
    {"mn",          DiskType::mn},
    {"exponential", DiskType::exponential}
  };

DiskType     dtype;
double       AA, HH;
double       ASHIFT;

// EXP support
//
#include <global.H>
#include <Particle.H>

double DiskDens(double R, double z, double phi)
{
  double ans = 0.0;
  static bool firsttime = true;

  switch (dtype) {
      

  case DiskType::constant:
    if (firsttime) std::cout << "Dens = constant" << std::endl;
    if (R < AA && fabs(z) < HH)
      ans = 1.0/(2.0*HH*M_PI*AA*AA);
    break;
  
  case DiskType::gaussian:
    if (firsttime) std::cout << "Dens = gaussian" << std::endl;
    if (fabs(z) < HH)
      ans = 1.0/(2.0*HH*2.0*M_PI*AA*AA)*
	exp(-R*R/(2.0*AA*AA));
    break;
    
  case DiskType::mn:
    if (firsttime) std::cout << "Dens = mn" << std::endl;
    {
      double Z2 = z*z + HH*HH;
      double Z  = sqrt(Z2);
      double Q  = AA + Z;
      ans = 0.25*HH*HH/M_PI*(AA*R*R + (AA + 3.0*Z)*Q*Q)/( pow(R*R + Q*Q, 2.5) * Z * Z2 );
    }
    break;

  default:
  case DiskType::exponential:
    if (firsttime) std::cout << "Dens = exponential" << std::endl;
    {
      double f = cosh(z/HH);
      ans = exp(-R/AA)/(4.0*M_PI*AA*AA*HH*f*f);
    }
    break;
  }

  firsttime = false;

  return ans;
}

double dcond(double R, double z, double phi, int M)
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
  double x = R*cos(phiS) - ASHIFT*AA;
  double y = R*sin(phiS);

  return DiskDens(sqrt(x*x + y*y), z, atan2(y, x));
}

// Mapping functions for integration
//
double r_to_x(double r, double a)
{
  return r/sqrt(r*r + a*a);
}

double x_to_r(double x, double a)
{
  return x*a/sqrt(1.0 - x*x);
}

double drdx(double x, double a)
{
  return a*pow(1.0 - x*x, -1.5);
}

int 
main(int ac, char **av)
{
  //====================
  // Inialize MPI stuff
  //====================
  
  local_init_mpi(ac, av);

  //====================
  // Begin opt parsing
  //====================

  int          LMAX;
  int          NMAX;
  int          NUMR;
  int          SCMAP;
  double       RCYLMIN;
  double       RCYLMAX;
  double       ASCALE;
  double       HSCALE;
  double       RFACTOR;
  double       AEXP;
  double       HEXP;
  int          RNUM;
  int          PNUM;
  int          TNUM;
  int          VFLAG;
  bool         expcond;
  bool         LOGR;
  bool         LOGR2;
  bool         CHEBY;
  int          CMAPR;
  int          CMAPZ;
  int          CMTYPE;
  int          NCHEB;
  int          NDR;
  int          NDZ;
  int          NHR;
  int          NHT;
  int          MMAX;
  int          NUMX;
  int          NUMY;
  int          NOUT;
  int          NFRC;
  int          NORDER;
  int          NODD;
  double       scale_height;
  double       scale_length;
  double       ppower;
  bool         SVD;
  int          NINT;
  bool         DENS;
  bool         ignore;
  bool         orthotst;
  string       cachefile;
  string       config;
  string       disktype;
  string       dmodel;
  string       mtype;
  
  cxxopts::Options options("testcoefs", "Check orthgonality for cylindrical basis");

  options.add_options()
    ("h,help", "Print this help message")
    ("T,template", "Write template options file with current and all default values")
    ("c,config", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("condition", "Condition EmpCylSL deprojection from specified disk model (EXP or MN)",
     cxxopts::value<string>(dmodel))
    ("LMAX", "Maximum harmonic order",
     cxxopts::value<int>(LMAX)->default_value("4"))
    ("NMAX", "Maximum radial order",
     cxxopts::value<int>(NMAX)->default_value("8"))
    ("MMAX", "Maximum azimuthal order",
     cxxopts::value<int>(MMAX)->default_value("y"))
    ("NUMX", "Size of the (mapped) cylindrical radial grid",
     cxxopts::value<int>(NUMX)->default_value("256"))
    ("NUMY", "Size of the (mapped) cylindrical vertical grid",
     cxxopts::value<int>(NUMY)->default_value("128"))
    ("NOUT", "Maximum radial order for diagnostic basis dump",
     cxxopts::value<int>(NOUT)->default_value("18")) 
    ("NINT", "Number of Gauss-Legendre knots",
     cxxopts::value<int>(NINT)->default_value("100"))
    ("NUMR", "Size of radial grid for Spherical SL",
     cxxopts::value<int>(NUMR)->default_value("2000"))
    ("NODD", "Number of vertically odd basis functions per harmonic order",
     cxxopts::value<int>(NODD)->default_value("6"))
    ("NORDER", "Total number of basis functions per harmonic order",
     cxxopts::value<int>(NORDER)->default_value("32"))
    ("NFRC", "Number of force evaluations in radius",
     cxxopts::value<int>(NFRC)->default_value("32"))
    ("RCYLMIN", "Minimum disk radius",
     cxxopts::value<double>(RCYLMIN)->default_value("0.001"))
    ("RCYLMAX", "Maximum disk radius",
     cxxopts::value<double>(RCYLMAX)->default_value("20.0"))
    ("ASCALE", "Radial scale length for disk basis construction",
     cxxopts::value<double>(ASCALE)->default_value("1.0"))
    ("ASHIFT", "Fraction of scale length for shift in conditioning function",
     cxxopts::value<double>(ASHIFT)->default_value("0.0"))
    ("HSCALE", "Vertical scale length for disk basis construction",
     cxxopts::value<double>(HSCALE)->default_value("0.1"))
    ("RFACTOR", "Disk radial scaling for deprojection computation",
     cxxopts::value<double>(RFACTOR)->default_value("1.0"))
    ("AEXP", "Radial scale length for disk basis test",
     cxxopts::value<double>(AEXP)->default_value("1.0"))
    ("HEXP", "Vertical scale length for disk basis test",
     cxxopts::value<double>(HEXP)->default_value("0.1"))
    ("RNUM", "Number of radial knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(RNUM)->default_value("200"))
    ("PNUM", "Number of azimthal knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(PNUM)->default_value("80"))
    ("TNUM", "Number of cos(theta) knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(TNUM)->default_value("80"))
    ("VFLAG", "Verbose reporting flags for EmpCylSL",
     cxxopts::value<int>(VFLAG)->default_value("0"))
    ("expcond", "EOF conditioning from analytic density profile",
     cxxopts::value<bool>(expcond)->default_value("true"))
    ("DENS", "Compute density basis functions",
     cxxopts::value<bool>(DENS)->default_value("true"))
    ("ignore", "Ignore parameters in cachefile and recompute basis",
     cxxopts::value<bool>(ignore)->default_value("false"))
    ("orthotst", "Compute orthgonality of basis functions",
     cxxopts::value<bool>(orthotst)->default_value("false"))
    ("LOGR", "Logarithmc mapping for radial profile in EmpCylSL",
     cxxopts::value<bool>(LOGR)->default_value("false"))
    ("LOGR2", "Logarithmc mapping for radial profile in coefficient evaluation",
     cxxopts::value<bool>(LOGR2)->default_value("false"))
    ("CMAPR", "Radial coordinate mapping type for cylindrical grid (0=none, 1=rational function)",
     cxxopts::value<int>(CMAPR)->default_value("1"))
    ("CMAPZ", "Vertical coordinate mapping type for cylindrical grid (0=none, 1=rational function)",
     cxxopts::value<int>(CMAPZ)->default_value("1"))
    ("scale_height", "Scale height for analytic disk density",
     cxxopts::value<double>(scale_height)->default_value("0.01"))
    ("scale_length", "Scale length for analytic disk density",
     cxxopts::value<double>(scale_length)->default_value("0.002"))
    ("PPOW", "Power-law exponent for analytic power law model",
     cxxopts::value<double>(ppower)->default_value("3.0"))
    ("cachefile", "Cache file for EmpCylSL grid",
     cxxopts::value<std::string>(cachefile)->default_value(".eof.cache.file"))
    ("disktype", "Desired disktype for analytic computation of basis (constant, gaussian, mn, exponential)",
     cxxopts::value<std::string>(disktype)->default_value("exponential"))
    ("mtype", "Desired deprojection sphericla models (Exponential, Gaussian, Plummer)",
     cxxopts::value<std::string>(mtype)->default_value("Exponential"))
     ;
  
  // Parse command line for control and critical parameters
  //
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << options.help() << std::endl << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  // Write YAML template config file and exit
  //
  if (vm.count("template")) {

    NOUT = std::min<int>(NOUT, NORDER);

    // Write template file
    //
    if (myid==0) SaveConfig(vm, options, "template.yaml");

    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the YAML config file
  //
  if (vm.count("input")) {
    try {
      vm = LoadConfig(options, config);
    } catch (cxxopts::OptionException& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }

  // Set EmpCylSL mtype
  //
  EmpCylSL::mtype = EmpCylSL::Exponential;
  if (vm.count("mtype")) {
    if (mtype.compare("Exponential")==0)
      EmpCylSL::mtype = EmpCylSL::Exponential;
    else if (mtype.compare("Gaussian")==0)
      EmpCylSL::mtype = EmpCylSL::Gaussian;
    else if (mtype.compare("Plummer")==0) {
      EmpCylSL::mtype = EmpCylSL::Plummer;
      EmpCylSL::PPOW  = ppower;
    } else {
      if (myid==0) std::cout << "No EmpCylSL EmpModel named <"
			     << mtype << ">, valid types are: "
			     << "Exponential, Gaussian, Plummer" << std::endl;
      MPI_Finalize();
      return -1;
    }
  }

  // Set DiskType
  //
  std::transform(disktype.begin(), disktype.end(), disktype.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  auto dit = dtlookup.find(disktype);

  if (dit == dtlookup.end()) {
    if (myid==0) {
      std::cout << "DiskType error in configuration file: "
		<< disktype << std::endl;
      std::cout << "Valid options are: ";
      for (auto v : dtlookup) std::cout << v.first << " ";
      std::cout << std::endl;
    }
    MPI_Finalize();
    return -1;
  }

  dtype = dit->second;

  // Report dtype
  //
  if (myid==0)
    std::cout << "DiskType is <" << disktype << ">" << std::endl;

  //====================
  // OpenMP control
  //====================

#ifdef HAVE_OMP_H
  omp_set_num_threads(nthrds);
#pragma omp parallel
  {
    int numthrd = omp_get_num_threads();
    int myid = omp_get_thread_num();
    if (myid==0)
      std::cout << "Number of threads=" << numthrd << std::endl;
  }
#endif

  //====================
  // Okay, now begin ...
  //====================

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  // set_fpu_handler();         // Make gdb trap FPU exceptions
  set_fpu_gdb_handler();	// Make gdb trap FPU exceptions
#endif
  
  //===========================Cylindrical expansion===========================


  // Set parameters from the given CACHEFILE
  //
  if (not ignore) {

    std::ifstream in(cachefile);
    if (!in) {
      std::cerr << "Error opening cachefile named <" 
		<< cachefile << "> . . ."
		<< std::endl
		<< "I will build <" << cachefile
		<< "> but it will take some time."
		<< std::endl
		<< "If this is NOT what you want, "
		<< "stop this routine and specify the correct file."
		<< std::endl;
    } else {

      // Attempt to read magic number
      //
      const unsigned int hmagic = 0xc0a57a1; // Basis magic #
      unsigned int tmagic;
      in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

      if (tmagic == hmagic) {

	// YAML size
	//
	unsigned ssize;
	in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

	// Make and read char buffer
	//
	auto buf = std::make_unique<char[]>(ssize+1);
	in.read(buf.get(), ssize);
	buf[ssize] = 0;		// Null terminate

	YAML::Node node;
      
	try {
	  node = YAML::Load(buf.get());
	}
	catch (YAML::Exception& error) {
	  if (myid==0)
	    std::cerr << "YAML: error parsing <" << buf.get() << "> "
		      << "in " << __FILE__ << ":" << __LINE__ << std::endl
		      << "YAML error: " << error.what() << std::endl;
	  throw error;
	}
	
	// Get parameters
	//
	MMAX   = node["mmax"  ].as<int>();
	NUMX   = node["numx"  ].as<int>();
	NUMY   = node["numy"  ].as<int>();
	NMAX   = node["nmax"  ].as<int>();
	NORDER = node["norder"].as<int>();
	DENS   = node["dens"  ].as<bool>();
	// RMIN   = node["rmin"  ].as<double>();
	// RMAX   = node["rmax"  ].as<double>();
	ASCALE = node["ascl"  ].as<double>();
	HSCALE = node["hscl"  ].as<double>();

	if (node["cmap"])		// Backwards compatibility
	  CMAPR  = node["cmap"  ].as<int>();
	else
	  CMAPR  = node["cmapr" ].as<int>();
	
	if (node["cmapz"])	// Backwards compatibility
	  CMAPZ  = node["cmapz" ].as<int>();

      } else {

	// Rewind file
	//
	in.clear();
	in.seekg(0);

      int tmp;
      
      in.read((char *)&MMAX,    sizeof(int));
      in.read((char *)&NUMX,    sizeof(int));
      in.read((char *)&NUMY,    sizeof(int));
      in.read((char *)&NMAX,    sizeof(int));
      in.read((char *)&NORDER,  sizeof(int));
      
      in.read((char *)&tmp,     sizeof(int)); 
      if (tmp) DENS = true;
      else     DENS = false;
      
      in.read((char *)&CMTYPE,  sizeof(int)); 
      in.read((char *)&RCYLMIN, sizeof(double));
      in.read((char *)&RCYLMAX, sizeof(double));
      in.read((char *)&ASCALE,  sizeof(double));
      in.read((char *)&HSCALE,  sizeof(double));
      }
    }
  }
    
  // Limit value of NOUT
  //
  NOUT = std::min<int>(NOUT, NORDER);


  // Assign parameters
  //
  EmpCylSL::RMIN        = RCYLMIN;
  EmpCylSL::RMAX        = RCYLMAX;
  EmpCylSL::NUMX        = NUMX;
  EmpCylSL::NUMY        = NUMY;
  EmpCylSL::NUMR        = NUMR;
  EmpCylSL::NOUT        = NOUT;
  EmpCylSL::CMAPR       = CMAPR;
  EmpCylSL::CMAPZ       = CMAPZ;
  EmpCylSL::VFLAG       = VFLAG;
  EmpCylSL::logarithmic = LOGR;
  EmpCylSL::DENS        = DENS;

				// For DiskDens
  AA = ASCALE;
  HH = HSCALE;

                                // Create expansion only if needed . . .
  bool save_eof = false;

  std::shared_ptr<EmpCylSL> expandd =
    std::make_shared<EmpCylSL>(NMAX, LMAX, MMAX, NORDER, ASCALE, HSCALE, NODD, cachefile);

  // Try to read existing cache to get EOF
  //
  if (not ignore) {
    save_eof = expandd->read_cache();
    if (myid==0) {
      if (save_eof) {
	std::cout << "EmpCylSL requested cache read: GOOD, continuing" << std::endl;
      } else {
	std::cout << "EmpCylSL requested cache read FAIL: exiting ..." << std::endl;
      }
    }
    if (not save_eof) {
      MPI_Finalize();
      return 4;
    }

  } else {

    // Use these user models to deproject for the EOF spherical basis
    //
    if (vm.count("condition")) {
      // The scale in EmpCylSL is assumed to be 1 so we compute the
      // height relative to the length
      //
      double H = scale_height/scale_length;
      
      // The model instance (you can add others in DiskModels.H
      //
      EmpCylSL::AxiDiskPtr model;
      
      if (dmodel.compare("MN")==0) // Miyamoto-Nagai
	model = std::make_shared<MNdisk>(1.0, H);
      else			// Default to exponential
	model = std::make_shared<Exponential>(1.0, H);
      
      expandd->create_deprojection(H, RFACTOR, NUMR, RNUM, model);
    }
    
    // Regenerate EOF from analytic density
    //
    if (expcond and not save_eof) {
      expandd->generate_eof(RNUM, PNUM, TNUM, dcond);
      save_eof = true;
    }
  }
  
  //===========================================================================
  // Compute coefficients
  //===========================================================================

  std::vector<double> coefs(NOUT, 0.0);
  LegeQuad lq(NINT);

  // Reassign DiskDens scales
  //
  AA = AEXP;
  HH = HEXP;

  std::cout << "A    = " << ASCALE  << std::endl
	    << "H    = " << HSCALE  << std::endl
	    << "AA   = " << AA      << std::endl
	    << "HH   = " << HH      << std::endl
	    << "Rmin = " << RCYLMIN << std::endl
	    << "Rmax = " << RCYLMAX << std::endl
	    << "Nord = " << NOUT    << std::endl;

  double xmin = r_to_x(RCYLMIN*AA, AA);
  double xmax = r_to_x(RCYLMAX*AA, AA);
  double ymax = r_to_x(RCYLMAX*AA, HH);


  std::map<std::pair<int, int>, double> orthochk;

  double totM = 0.0;

  if (LOGR2) {
    
    double Rmin = log(RCYLMIN*AA);
    double Rmax = log(RCYLMAX*AA);
    double Zmin = log(RCYLMIN*HH);
    double Zmax = log(RCYLMAX*AA);

    for (int i=0; i<NINT; i++) {	// Radial

      double x = Rmin + (Rmax - Rmin) * lq.knot(i);
      double R = exp(x);

      double facX = lq.weight(i) * 2.0*M_PI * R * R * (Rmax - Rmin);

      for (int j=0; j<NINT; j++) { // Vertical
	
	double y = Zmin + (Zmax - Zmin) * lq.knot(j);
	double z = exp(y);
	
	double fac = facX * lq.weight(j) * z * (Zmax - Zmin);
	double den = DiskDens(R, z, 0.0);
	totM += 2.0 * fac * den;

	fac *= -1.0;

	for (int n=0; n<NOUT; n++) {
	  double p, p2, d, d2, fr, fz, fp;
	  expandd->get_all(0, n, R, z, 0.0, p, d, fr, fz, fp);
	  coefs[n] += fac * p * den * 4.0*M_PI;
	  
	  if (orthotst) {
	    for (int n2=n; n2<NOUT; n2++) {
	      if (n2>n) expandd->get_all(0, n2, R, z, 0.0, p2, d2, fr, fz, fp);
	      else      d2 = d;
	      orthochk[{n, n2}] += fac * p * d2 * 4.0*M_PI;
	    }
	  }
	}
	
	for (int n=0; n<NOUT; n++) {
	  double p, p2, d, d2, fr, fz, fp;
	  expandd->get_all(0, n, R, -z, 0.0, p, d, fr, fz, fp);
	  coefs[n] += fac * p * den * 4.0*M_PI;

	  if (orthotst) {
	    for (int n2=n; n2<NOUT; n2++) {
	      if (n2>n) expandd->get_all(0, n2, R, -z, 0.0, p2, d2, fr, fz, fp);
	      else      d2 = d;
	      orthochk[{n, n2}] += fac * p * d2 * 4.0*M_PI;
	    }
	  }
	}
      }
    }

  } else {

    for (int i=0; i<NINT; i++) {	// Radial

      double x = xmin + (xmax - xmin) * lq.knot(i);
      double R = x_to_r(x, AA);

      double facX = lq.weight(i) * 2.0 * M_PI * R * drdx(x, AA) * (xmax - xmin);

      for (int j=0; j<NINT; j++) { // Vertical
	
	double y = ymax*(2.0*lq.knot(j) - 1.0);
	double z = x_to_r(y, HH);
	
	double fac = facX * lq.weight(j) * drdx(y, HH) * 2.0*ymax;

	double den = DiskDens(R, z, 0.0);
	totM += fac * den;

	fac *= -1.0;

	for (int n=0; n<NOUT; n++) {
	  double p, p2, d, d2, fr, fz, fp;
	  expandd->get_all(0, n, R, z, 0.0, p, d, fr, fz, fp);
	  coefs[n] += fac * p * den * 4.0*M_PI;

	  if (orthotst) {
	    for (int n2=n; n2<NOUT; n2++) {
	      if (n2>n) expandd->get_all(0, n2, R, z, 0.0, p2, d2, fr, fz, fp);
	      else      d2 = d;
	      orthochk[{n, n2}] += fac * p * d2 * 4.0*M_PI;
	    }
	  }
	}
      }
    }
  }

  std::cout << std::endl << "Total integrated mass: " << totM
	    << std::endl << std::endl;

  double cum = 0.0;
  for (int n=0; n<NOUT; n++) {
    cum += coefs[n]*coefs[n];
    std::cout << std::setw( 8) << n
	      << std::setw(18) << coefs[n]
	      << std::setw(18) << coefs[n]*coefs[n]
	      << std::setw(18) << cum
	      << std::endl;
  }

  std::cout << std::endl;

  for (auto v : orthochk)
    std::cout << std::setw( 3) << v.first.first
	      << std::setw( 3) << v.first.second
	      << std::setw(18) << v.second
	      << std::endl;

  // Set coefficients from std::vectors
  //
  std::vector<double> zero(NOUT, 0.0);
  expandd->set_coefs(0, coefs, zero, true);


  // Quick radial force check
  //
  double dx   = (xmax - xmin)/(NFRC-1);
  double z    = 0.0;
  double phi  = 0.0;
  double mass = 1.0;
  
  std::ofstream fout("testcoefs.compare");

  int nmin = std::min<int>(NOUT, 5);

  for (int j=0; j<NFRC; j++) {
    std::vector<double> dd(nmin);
    double p0, p, fr, fz, fp, d;
    double r = x_to_r(xmin + dx*j, AA);

    expandd->accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
    expandd->accumulated_dens_eval(r, z, phi, d);

    // Get density for n=0, 1, ... , nmin
    {
      double p1, fr1, fz1, fp1;	// Dummy variables
      for (int nn=0; nn<nmin; nn++) 
	expandd->get_all(0, nn, r, z, 0.0, p1, dd[nn], fr1, fz1, fp1);
    }

    
    double D, P, FR;

    if (dmodel.compare("MN")==0) { // Miyamoto-Nagai
      double zb = sqrt( z*z + HH*HH );
      double ab = AA + zb;
      double dn = sqrt( r*r + ab*ab );
      
      P  = -mass/dn;
      // double FZ = -mass*z*ab/(zb*dn*dn*dn);
      FR = -mass*r/(dn*dn*dn);
      D  = 0.25*HH*HH/M_PI*(AA*r*r + (AA + 3.0*zb)*ab*ab)/( pow(r*r + ab*ab, 2.5) * zb*zb*zb );

    } else {			// Default to exponential
      double y = r/(2.0*AA);
      double i0 = std::cyl_bessel_i(0, y);
      double k0 = std::cyl_bessel_k(0, y);
      double i1 = std::cyl_bessel_i(1, y);
      double k1 = std::cyl_bessel_k(1, y);
      P  = -0.5*mass*r/(AA*AA) * (i0*k1 - i1*k0);
      FR = -y/(AA*AA)*(i0*k0 - i1*k1);

      double f = cosh(z/HH);
      D  = exp(-r/AA)/(4.0*M_PI*AA*AA*HH*f*f);
    }

    fout << std::setw(18) << r	          // 1
	 << std::setw(18) << fr	          // 2
	 << std::setw(18) << FR	          // 3
	 << std::setw(18) << (fr - FR)/FR // 4
	 << std::setw(18) << p		  // 5
	 << std::setw(18) << P		  // 6
	 << std::setw(18) << (p - P)/P	  // 7
	 << std::setw(18) << d		  // 8
	 << std::setw(18) << D;		  // 9
    for (int nn=0; nn<nmin; nn++)
      fout << std::setw(18) << dd[nn];	  // 10+nn
    fout << std::endl;
  }


  //===========================================================================
  // shutdown MPI
  //===========================================================================

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

