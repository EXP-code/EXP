/*
  Check coefficients for cylindrical basis
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
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

// MDW classes
//
#include <numerical.h>
#include <gaussQ.h>
#include <isothermal.h>
#include <hernquist.h>
#include <model3d.h>
#include <biorth.h>
#include <SphericalSL.h>
#include <interp.h>
#include <EmpCylSL.h>
#include <DiskModels.H>

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
#include "SphericalSL.h"
#include "DiskHalo2.h" 
#include "localmpi.h"


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
double       ASCALE;
double       ASHIFT;
double       HSCALE;

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

  switch (dtype) {
      
  constant:
    if (R < ASCALE && fabs(z) < HSCALE)
      ans = 1.0/(2.0*HSCALE*M_PI*ASCALE*ASCALE);
    break;
  
  gaussian:
    if (fabs(z) < HSCALE)
      ans = 1.0/(2.0*HSCALE*2.0*M_PI*ASCALE*ASCALE)*
	exp(-R*R/(2.0*ASCALE*ASCALE));
    break;
    
  mn:
    {
      double Z2 = z*z + HSCALE*HSCALE;
      double Z  = sqrt(Z2);
      double Q  = ASCALE + Z;
      ans = 0.25*HSCALE*HSCALE/M_PI*(ASCALE*R*R + (ASCALE + 3.0*Z)*Q*Q)/( pow(R*R + Q*Q, 2.5) * Z * Z2 );
    }
    break;

  default:
  exponential:
    {
      double f = cosh(z/HSCALE);
      ans = exp(-R/ASCALE)/(4.0*M_PI*ASCALE*ASCALE*HSCALE*f*f);
    }
    break;
  }

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
  double x = R*cos(phiS) - ASHIFT*ASCALE;
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
  double       ZMIN;
  double       ZMAX;
  double       RMIN;
  double       RMAX;
  double       RCYLMIN;
  double       RCYLMAX;
  int          RNUM;
  int          PNUM;
  int          TNUM;
  int          VFLAG;
  bool         expcond;
  bool         LOGR;
  bool         CHEBY;
  int          CMAP;
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
  double       scale_height;
  double       scale_length;
  bool         SVD;
  int          NINT;
  bool         DENS;
  bool         ignore;
  string       cachefile;
  string       config;
  string       dtype;
  string       dmodel;
  string       mtype;
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("conf,c",          po::value<string>(&config),                                     "Write template options file with current and all default values")
    ("input,f",         po::value<string>(&config),                                     "Parameter configuration file")
    ("condition",       po::value<string>(&dmodel),                                     "Condition EmpCylSL deprojection from specified disk model (EXP or MN)")
    ("NINT",            po::value<int>(&NINT)->default_value(40),                       "Number of Gauss-Legendre knots")
    ("NUMR",            po::value<int>(&NUMR)->default_value(2000),                     "Size of radial grid for Spherical SL")
    ("RMIN",            po::value<double>(&RMIN)->default_value(0.0),                   "Minimum cylindrical radius")
    ("RMAX",            po::value<double>(&RMAX)->default_value(0.1),                   "Maximum cylindrical radius")
    ("ZMIN",            po::value<double>(&ZMIN)->default_value(0.0),                   "Minimum vertical height")
    ("ZMAX",            po::value<double>(&ZMAX)->default_value(0.1),                   "Maximum vertical height")
    ("RCYLMIN",         po::value<double>(&RCYLMIN)->default_value(0.001),              "Minimum disk radius")
    ("RCYLMAX",         po::value<double>(&RCYLMAX)->default_value(20.0),               "Maximum disk radius")
    ("ASCALE",          po::value<double>(&ASCALE)->default_value(1.0),                 "Radial scale length for disk basis construction")
    ("ASHIFT",          po::value<double>(&ASHIFT)->default_value(0.0),                 "Fraction of scale length for shift in conditioning function")
    ("HSCALE",          po::value<double>(&HSCALE)->default_value(0.1),                 "Vertical scale length for disk basis construction")
    ("RNUM",            po::value<int>(&RNUM)->default_value(200),                      "Number of radial knots for EmpCylSL basis construction quadrature")
    ("PNUM",            po::value<int>(&PNUM)->default_value(80),                       "Number of azimthal knots for EmpCylSL basis construction quadrature")
    ("TNUM",            po::value<int>(&TNUM)->default_value(80),                       "Number of cos(theta) knots for EmpCylSL basis construction quadrature")
    ("CMAP",            po::value<int>(&CMAP)->default_value(2),                        "Map coordinates from radius to tabled grid")
    ("SVD",             po::value<bool>(&SVD)->default_value(false),                    "Use svd for symmetric eigenvalue problesm")
    ("LOGR",            po::value<bool>(&LOGR)->default_value(false),                   "Make a logarithmic coordinate mapping")
    ("LMAX",            po::value<int>(&LMAX)->default_value(6),                        "Number of harmonics for Spherical SL for halo/spheroid")
    ("NMAX",            po::value<int>(&NMAX)->default_value(12),                       "Number of radial basis functions in Spherical SL for halo/spheroid")
    ("MMAX",            po::value<int>(&MMAX)->default_value(4),                        "Number of azimuthal harmonics for disk basis")
    ("NUMX",            po::value<int>(&NUMX)->default_value(256),                      "Radial grid size for disk basis table")
    ("NUMY",            po::value<int>(&NUMY)->default_value(128),                      "Vertical grid size for disk basis table")
    ("NORDER",          po::value<int>(&NORDER)->default_value(16),                     "Number of disk basis functions per M-order")
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
    ("mtype",           po::value<string>(&mtype),                                              "Spherical deprojection model for EmpCylSL (one of: Exponential, Gaussian, Plummer)")
    ("DTYPE",           po::value<string>(&dtype)->default_value("exponential"),                 "Disk type for condition (one of: constant, gaussian, mn, exponential)")
    ("ignore",          po::value<bool>(&ignore)->default_value(false),                 "Ignore any existing cache file and recompute the EOF")
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

    NOUT = std::min<int>(NOUT, NORDER);

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
    else if (mtype.compare("Plummer")==0)
      EmpCylSL::mtype = EmpCylSL::Plummer;
    else {
      if (myid==0) std::cout << "No EmpCylSL EmpModel named <"
			     << mtype << ">, valid types are: "
			     << "Exponential, Gaussian, Plummer" << std::endl;
      MPI_Finalize();
      return -1;
    }
  }

  // Set DiskType
  //
  std::transform(dtype.begin(), dtype.end(), dtype.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  auto dit = dtlookup.find(dtype);

  if (dit == dtlookup.end()) {
    if (myid==0) {
      std::cout << "DiskType error in configuration file: "
		<< dtype << std::endl;
      std::cout << "Valid options are: ";
      for (auto v : dtlookup) std::cout << v.first << " ";
      std::cout << std::endl;
    }
    MPI_Finalize();
    return -1;
  }

  dtype = dit->second;

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

      int tmp;
    
      in.read((char *)&MMAX,    sizeof(int));
      in.read((char *)&NUMX,    sizeof(int));
      in.read((char *)&NUMY,    sizeof(int));
      in.read((char *)&NMAX,    sizeof(int));
      in.read((char *)&NORDER,  sizeof(int));
      
      in.read((char *)&tmp,     sizeof(int)); 
      if (tmp) DENS = true;
      else     DENS = false;
      
      in.read((char *)&CMAP,    sizeof(int)); 
      in.read((char *)&RCYLMIN, sizeof(double));
      in.read((char *)&RCYLMAX, sizeof(double));
      in.read((char *)&ASCALE,  sizeof(double));
      in.read((char *)&HSCALE,  sizeof(double));
    }
  }

  EmpCylSL::RMIN        = RCYLMIN;
  EmpCylSL::RMAX        = RCYLMAX;
  EmpCylSL::NUMX        = NUMX;
  EmpCylSL::NUMY        = NUMY;
  EmpCylSL::NUMR        = NUMR;
  EmpCylSL::NOUT        = NOUT;
  EmpCylSL::CMAP        = CMAP;
  EmpCylSL::VFLAG       = VFLAG;
  EmpCylSL::logarithmic = LOGR;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::USESVD      = SVD;
  EmpCylSL::CACHEFILE   = cachefile;

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
      
      expandd->create_deprojection(H, NUMR, RNUM, model);
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

  std::vector<double> coefs(NORDER, 0.0);
  LegeQuad lq(NINT);

  std::cout << "A   =" << ASCALE << std::endl
	    << "H   =" << HSCALE << std::endl
	    << "Rmin=" << RCYLMIN << std::endl
	    << "Rmax=" << RCYLMAX << std::endl
	    << "RMIN=" << RMIN << std::endl
	    << "RMAX=" << RMAX << std::endl;

  double xmin = r_to_x(RMIN, ASCALE);
  double xmax = r_to_x(RMAX, ASCALE);
  double ymin = r_to_x(ZMIN, HSCALE);
  double ymax = r_to_x(ZMAX, HSCALE);

  for (int i=1; i<=NINT; i++) {	// Radial

    double x = xmin + (xmax - xmin) * lq.knot(i);
    double R = x_to_r(x, ASCALE);

    double facX = lq.weight(i) * 2.0 * M_PI * R * drdx(x, ASCALE) * (xmax - xmin);

    for (int j=1; j<=NINT; j++) { // Vertical

      double y = ymin + (ymax - ymin)*lq.knot(j);
      double z = x_to_r(y, HSCALE);

      double fac = facX * lq.weight(j) * 2.0 * drdx(y, HSCALE) * (ymax - ymin);

      for (int n=0; n<NORDER; n++) {
	double p, d, fr, fz, fp;
	expandd->get_all(0, n, R, z, 0.0, p, d, fr, fz, fp);
	coefs[n] += fac * p * DiskDens(R, z, 0.0) * 4.0*M_PI;
      }
    }
  }

  double cum = 0.0;
  for (int n=0; n<NORDER; n++) {
    cum += coefs[n]*coefs[n];
    std::cout << std::setw( 8) << n
	      << std::setw(18) << coefs[n]
	      << std::setw(18) << coefs[n]*coefs[n]
	      << std::setw(18) << cum
	      << std::endl;
  }

  // Set coefficients from std::vectors
  //
  std::vector<double> zero(NORDER, 0.0);
  expandd->set_coefs(0, coefs, zero, true);


  // Quick radial force check
  //
  double dr   = (RMAX - RMIN)/(NFRC-1);
  double z    = 0.0;
  double phi  = 0.0;
  double mass = 1.0;
  
  std::ofstream fout("testcoefs.compare");

  for (int j=0; j<NFRC; j++) {
    double p0, p, fr, fz, fp;
    double r = RMIN + dr*j;

    expandd->accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
    
    double zb = sqrt( z*z + HSCALE*HSCALE );
    double ab = ASCALE + zb;
    double dn = sqrt( r*r + ab*ab );
    
    double PP = -mass/dn;
    double FR = -mass*r/(dn*dn*dn);
    double FZ = -mass*z*ab/(zb*dn*dn*dn);

    fout << std::setw(18) << r
	 << std::setw(18) << -fr
	 << std::setw(18) << FR
	 << std::setw(18) << (-fr - FR)/FR
	 << std::setw(18) << p
	 << std::setw(18) << 1.0/dn
	 << std::setw(18) << (p - 1.0/dn)*dn
	 << std::endl;
  }


  //===========================================================================
  // shutdown MPI
  //===========================================================================

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

