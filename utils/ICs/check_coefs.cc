/*
  Check orthgonality for cylindrical basis
*/
                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <fenv.h>

// Boost stuff
//
#include <boost/math/special_functions/bessel.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>
#include <boost/make_unique.hpp>

#include <yaml-cpp/yaml.h>	// YAML support

namespace po = boost::program_options;

#include <config.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

// MDW classes
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

#include <Particle.H>

int VERBOSE        = 4;
int nthrds         = 1;
int this_step      = 0;
unsigned multistep = 0;
unsigned maxlev    = 100;
int mstep          = 1;
int Mstep          = 1;
char threading_on  = 0;
double tpos        = 0.0;
double tnow        = 0.0;

vector<int> stepL(1, 0), stepN(1, 1);
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
string outdir, runtag;
  
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
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("conf,c",          po::value<string>(&config),                                     "Write template options file with current and all default values")
    ("input,f",         po::value<string>(&config),                                     "Parameter configuration file")
    ("condition",       po::value<string>(&dmodel),                                     "Condition EmpCylSL deprojection from specified disk model (EXP or MN)")
    ("NINT",            po::value<int>(&NINT)->default_value(100),                       "Number of Gauss-Legendre knots")
    ("NUMR",            po::value<int>(&NUMR)->default_value(2000),                     "Size of radial grid for Spherical SL")
    ("RCYLMIN",         po::value<double>(&RCYLMIN)->default_value(0.001),              "Minimum disk radius")
    ("RCYLMAX",         po::value<double>(&RCYLMAX)->default_value(20.0),               "Maximum disk radius")
    ("ASCALE",          po::value<double>(&ASCALE)->default_value(1.0),                 "Radial scale length for disk basis construction")
    ("ASHIFT",          po::value<double>(&ASHIFT)->default_value(0.0),                 "Fraction of scale length for shift in conditioning function")
    ("HSCALE",          po::value<double>(&HSCALE)->default_value(0.1),                 "Vertical scale length for disk basis construction")
    ("RFACTOR",         po::value<double>(&RFACTOR)->default_value(1.0),                "Disk radial scaling for deprojection computation")
    ("AEXP",            po::value<double>(&AEXP)->default_value(1.0),                 "Radial scale length for disk basis test")
    ("HEXP",            po::value<double>(&HEXP)->default_value(0.1),                 "Vertical scale length for disk basis test")
    ("RNUM",            po::value<int>(&RNUM)->default_value(200),                      "Number of radial knots for EmpCylSL basis construction quadrature")
    ("PNUM",            po::value<int>(&PNUM)->default_value(80),                       "Number of azimthal knots for EmpCylSL basis construction quadrature")
    ("TNUM",            po::value<int>(&TNUM)->default_value(80),                       "Number of cos(theta) knots for EmpCylSL basis construction quadrature")
    ("CMAPR",           po::value<int>(&CMAPR)->default_value(1),                     "Radial coordinate mapping type for cylindrical grid (0=none, 1=rational fct)")
    ("CMAPZ",           po::value<int>(&CMAPZ)->default_value(1),                     "Vertical coordinate mapping type for cylindrical grid (0=none, 1=sech, 2=power in z")
    ("SVD",             po::value<bool>(&SVD)->default_value(false),                    "Use svd for symmetric eigenvalue problesm")
    ("LOGR",            po::value<bool>(&LOGR)->default_value(false),                   "Make a logarithmic coordinate mapping")
    ("LOGR2",           po::value<bool>(&LOGR2)->default_value(false),                   "Make a logarithmic coordinate mapping for coefficient eval")
    ("LMAX",            po::value<int>(&LMAX)->default_value(6),                        "Number of harmonics for Spherical SL for halo/spheroid")
    ("NMAX",            po::value<int>(&NMAX)->default_value(12),                       "Number of radial basis functions in Spherical SL for halo/spheroid")
    ("MMAX",            po::value<int>(&MMAX)->default_value(4),                        "Number of azimuthal harmonics for disk basis")
    ("NUMX",            po::value<int>(&NUMX)->default_value(256),                      "Radial grid size for disk basis table")
    ("NUMY",            po::value<int>(&NUMY)->default_value(128),                      "Vertical grid size for disk basis table")
    ("NORDER",          po::value<int>(&NORDER)->default_value(16),                     "Number of disk basis functions per M-order")
    ("NODD",            po::value<int>(&NODD)->default_value(-1),                       "Number of antisymmetric vertical functions per M-order")
    ("NOUT",            po::value<int>(&NOUT)->default_value(1000),                     "Number of radial basis functions to output for each harmonic order")
    ("NFRC",            po::value<int>(&NFRC)->default_value(100),                     "Number of radial knots for force check")
    ("DENS",            po::value<bool>(&DENS)->default_value(true),                    "Compute the density basis functions")
    ("VFLAG",           po::value<int>(&VFLAG)->default_value(0),                       "Output flags for EmpCylSL")
    ("threads",         po::value<int>(&nthrds)->default_value(1),                      "Number of lightweight threads")
    ("expcond",         po::value<bool>(&expcond)->default_value(true),                 "Use analytic density function for computing EmpCylSL basis")
    ("cachefile",       po::value<string>(&cachefile)->default_value(".eof.cache.file"),        "Name of EOF cache file")
    ("runtag",          po::value<string>(&runtag)->default_value("run000"),                    "Label prefix for diagnostic images")
    ("scale_height",    po::value<double>(&scale_height)->default_value(0.1),           "Scale height for disk realization")
    ("scale_length",    po::value<double>(&scale_length)->default_value(2.0),           "Scale length for disk realization")
    ("mtype",           po::value<string>(&mtype),                                              "Spherical deprojection model for EmpCylSL (one of: Exponential, Gaussian, Plummer, Power)")
    ("DTYPE",           po::value<string>(&disktype)->default_value("exponential"),             "Disk type for condition (one of: constant, gaussian, mn, exponential)")
    ("PPOW",            po::value<double>(&ppower)->default_value(5.0),             "Power-law density exponent for general Plummer density for EMP construction")
    ("ignore",          po::value<bool>(&ignore)->default_value(false),                 "Ignore any existing cache file and recompute the EOF")
    ("ortho",           po::value<bool>(&orthotst)->default_value(false),               "Check basis orthogonality by scalar product")
    ;
        
  po::variables_map vm;
  
  // Parse command line for control and critical parameters
  //
  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error on command line: "
			   << e.what() << std::endl;
    MPI_Finalize();
    return -1;
  }
  
  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      const char *mesg = "Generates a Monte Carlo realization of a halo\nwith an embedded disk using Jeans' equations.";
      std::cout << mesg << std::endl
		<< desc << std::endl << std::endl
		<< "Examples: " << std::endl
		<< "\t" << "Use parameters read from a config file in INI style"  << std::endl
		<< "\t" << av[0] << " --input=gendisk.config"  << std::endl << std::endl
		<< "\t" << "Generate a template config file in INI style from current defaults"  << std::endl
		<< "\t" << av[0] << " --conf=template.config" << std::endl << std::endl
		<< "\t" << "Override a single parameter in a config file from the command line"  << std::endl
		<< "\t" << av[0] << "--LMAX=8 --conf=template.config" << std::endl << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  // Write template config file in INI style and exit
  //
  if (vm.count("conf")) {
    // Do not overwrite existing config file
    //
    if (boost::filesystem::exists(config)) {
      if (myid == 0)
	std::cerr << av[0] << ": config file <" << config
		  << "> exists, will not overwrite" << std::endl;
      MPI_Finalize();
      return 2;
    }

    // Write template file
    //
    if (myid==0) {
      std::ofstream out(config);

      if (out) {
	// Iterate map and print out key--value pairs and description
	//
	for (const auto& it : vm) {
				// Don't write this parameter
	  if (it.first.find("conf")==0) continue;

	  out << std::setw(20) << std::left << it.first << " = ";
	  auto& value = it.second.value();
	  if (auto v = boost::any_cast<uint32_t>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<int>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<unsigned>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<float>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<double>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<bool>(&value))
	    out << std::setw(32) << std::left << std::boolalpha << *v;
	  else if (auto v = boost::any_cast<std::string>(&value))
	    out << std::setw(32) << std::left << *v;
	  else
	    out << "error";


	  //                               NO approximations -----+
	  // Add description as a comment                         |
	  //                                                      V
	  const po::option_description& rec = desc.find(it.first, false);
	  out << " # " << rec.description() << std::endl;
	}
      } else {
	if (myid==0)
	  std::cerr << av[0] << ": error opening template config file <"
		    << config << ">" << std::endl;
      }
    }
    MPI_Finalize();
    return 3;
  }

  // Read parameters fron the config file
  //
  if (vm.count("input")) {
    try {
      std::ifstream in(config);
      po::store(po::parse_config_file(in, desc), vm);
      po::notify(vm);    
    } catch (po::error& e) {
      if (myid==0) std::cout << "Option error in configuraton file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return -1;
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
	auto buf = boost::make_unique<char[]>(ssize+1);
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
  EmpCylSL::USESVD      = SVD;
  EmpCylSL::CACHEFILE   = cachefile;

				// For DiskDens
  AA = ASCALE;
  HH = HSCALE;

                                // Create expansion only if needed . . .
  bool save_eof = false;

  boost::shared_ptr<EmpCylSL> expandd =
    boost::make_shared<EmpCylSL>(NMAX, LMAX, MMAX, NORDER, ASCALE, HSCALE);

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
	model = boost::make_shared<MNdisk>(1.0, H);
      else			// Default to exponential
	model = boost::make_shared<Exponential>(1.0, H);
      
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

    for (int i=1; i<=NINT; i++) {	// Radial

      double x = Rmin + (Rmax - Rmin) * lq.knot(i);
      double R = exp(x);

      double facX = lq.weight(i) * 2.0*M_PI * R * R * (Rmax - Rmin);

      for (int j=1; j<=NINT; j++) { // Vertical
	
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

    for (int i=1; i<=NINT; i++) {	// Radial

      double x = xmin + (xmax - xmin) * lq.knot(i);
      double R = x_to_r(x, AA);

      double facX = lq.weight(i) * 2.0 * M_PI * R * drdx(x, AA) * (xmax - xmin);

      for (int j=1; j<=NINT; j++) { // Vertical
	
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
      double i0 = boost::math::cyl_bessel_i(0, y);
      double k0 = boost::math::cyl_bessel_k(0, y);
      double i1 = boost::math::cyl_bessel_i(1, y);
      double k1 = boost::math::cyl_bessel_k(1, y);
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

