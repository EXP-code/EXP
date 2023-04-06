/*
  Generates a cylindrical basis cache file
*/

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <random>
#include <vector>
#include <cmath>

#include <fenv.h>

#include <config_exp.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

// EXP classes
//
#include <libvars.H>
#include <numerical.H>
#include <gaussQ.H>
#include <isothermal.H>
#include <hernquist_model.H>
#include <model3d.H>
#include <biorth.H>
#include <SphericalSL.H>
#include <interp.H>
#include <EmpCylSL.H>
#include <DiskModels.H>
#include <cxxopts.H>		// Command-line parsing
#include <EXPini.H>		// Ini-style config

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


// Global variables
//
enum DiskType { constant, gaussian, mn, exponential, doubleexpon };

std::map<std::string, DiskType> dtlookup =
  { {"constant",    DiskType::constant},
    {"gaussian",    DiskType::gaussian},
    {"mn",          DiskType::mn},
    {"exponential", DiskType::exponential},
    {"doubleexpon", DiskType::doubleexpon}
  };

DiskType     DTYPE;
double       ASCALE;
double       ASHIFT;
double       HSCALE;
double       RTRUNC  = 1.0;
double       RWIDTH  = 0.0;
double       ARATIO  = 1.0;
double       HRATIO  = 1.0;
double       DWEIGHT = 1.0;

#include <Particle.H>

double DiskDens(double R, double z, double phi)
{
  double ans = 0.0;

  switch (DTYPE) {

  case DiskType::constant:
    if (R < ASCALE && fabs(z) < HSCALE)
      ans = 1.0/(2.0*HSCALE*M_PI*ASCALE*ASCALE);
    break;

  case DiskType::gaussian:
    if (fabs(z) < HSCALE)
      ans = 1.0/(2.0*HSCALE*2.0*M_PI*ASCALE*ASCALE)*
	exp(-R*R/(2.0*ASCALE*ASCALE));
    break;

  case DiskType::mn:
    {
      double Z2 = z*z + HSCALE*HSCALE;
      double Z  = sqrt(Z2);
      double Q2 = (ASCALE + Z)*(ASCALE + Z);
      ans = 0.25*HSCALE*HSCALE/M_PI*(ASCALE*R*R + (ASCALE + 3.0*Z)*Q2)/( pow(R*R + Q2, 2.5) * Z*Z2 );
    }
    break;

  case DiskType::doubleexpon:
    {
      double a1 = ASCALE;
      double a2 = ASCALE*ARATIO;
      double h1 = HSCALE;
      double h2 = HSCALE*HRATIO;
      double w1 = 1.0/(1.0+DWEIGHT);
      double w2 = DWEIGHT/(1.0+DWEIGHT);

      double f1 = cosh(z/h1);
      double f2 = cosh(z/h2);

      ans =
	w1*exp(-R/a1)/(4.0*M_PI*a1*a1*h1*f1*f1) +
	w2*exp(-R/a2)/(4.0*M_PI*a2*a2*h2*f2*f2) ;
    }
    break;
  case DiskType::exponential:
  default:
    {
      double f = cosh(z/HSCALE);
      ans = exp(-R/ASCALE)/(4.0*M_PI*ASCALE*ASCALE*HSCALE*f*f);
    }
    break;
  }

  if (RWIDTH>0.0) ans *= erf((RTRUNC-R)/RWIDTH);

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

  double       RCYLMIN;
  double       RCYLMAX;
  double       RFACTOR;
  int          RNUM;
  int          PNUM;
  int          TNUM;
  int          VFLAG;
  bool         expcond;
  bool         LOGR;
  int          NUMR;
  int          CMAPR;
  int          CMAPZ;
  int          CMTYPE;
  int          NMAX2;
  int          LMAX2;
  int          MMAX;
  int          NUMX;
  int          NUMY;
  int          NOUT;
  int          NORDER;
  int          NODD;
  double       PPower;
  bool         DENS;
  std::string  cachefile;
  std::string  config;
  std::string  dtype;
  std::string  dmodel;
  std::string  mtype;
  std::string  ctype;

  const std::string mesg("Generates an EmpCylSL cache\n");

  cxxopts::Options options(av[0], mesg);

  options.add_options()
    ("h,help", "Print this help message")
    ("T,template", "Write template options file with current and all default values")
    ("c,config", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("deproject", "The EmpCylSL deprojection from specified disk model (EXP or MN)",
     cxxopts::value<string>(dmodel)->default_value("EXP"))
    ("cachefile", "The cache file for the cylindrical basis",
     cxxopts::value<string>(cachefile)->default_value(".eof.cache.file"))
    ("ctype", "DiskHalo radial coordinate scaling type (one of: Linear, Log, Rat)",
     cxxopts::value<string>(ctype)->default_value("Log"))
    ("LMAX2", "Maximum angular order for spherical basis in adaptive construction of the cylindrical basis",
     cxxopts::value<int>(LMAX2)->default_value("48"))
    ("NMAX2", "Maximum radial order for the spherical basis in adapative construction of the cylindrical basis",
     cxxopts::value<int>(NMAX2)->default_value("48"))
    ("MMAX", "Maximum azimuthal order for the cylindrical basis",
     cxxopts::value<int>(MMAX)->default_value("6"))
    ("NUMX", "Size of the (mapped) cylindrical radial grid",
     cxxopts::value<int>(NUMX)->default_value("256"))
    ("NUMY", "Size of the (mapped) cylindrical vertical grid",
     cxxopts::value<int>(NUMY)->default_value("128"))
    ("dtype", "Spherical model type for adpative basis creation",
     cxxopts::value<string>(dtype)->default_value("exponential"))
    ("NOUT", "Maximum radial order for diagnostic basis dump",
     cxxopts::value<int>(NOUT)->default_value("18"))
    ("NODD", "Number of vertically odd basis functions per harmonic order",
     cxxopts::value<int>(NODD)->default_value("6"))
    ("NORDER", "Total number of basis functions per harmonic order",
     cxxopts::value<int>(NORDER)->default_value("32"))
    ("VFLAG", "Diagnostic flag for EmpCylSL",
     cxxopts::value<int>(VFLAG)->default_value("31"))
    ("expcond", "Use analytic target density rather than particle distribution",
     cxxopts::value<bool>(expcond)->default_value("true"))
    ("LOGR", "Logarithmic scaling for model table in EmpCylSL",
     cxxopts::value<bool>(LOGR)->default_value("true"))
    ("DENS", "Compute and cache basis density field",
     cxxopts::value<bool>(DENS)->default_value("true"))
    ("RCYLMIN", "Minimum disk radius for EmpCylSL",
     cxxopts::value<double>(RCYLMIN)->default_value("0.001"))
    ("RCYLMAX", "Maximum disk radius for EmpCylSL",
     cxxopts::value<double>(RCYLMAX)->default_value("20.0"))
    ("ASCALE", "Radial scale length for disk basis construction",
     cxxopts::value<double>(ASCALE)->default_value("1.0"))
    ("ASHIFT", "Fraction of scale length for shift in conditioning function",
     cxxopts::value<double>(ASHIFT)->default_value("0.0"))
    ("HSCALE", "Vertical scale length for disk basis construction",
     cxxopts::value<double>(HSCALE)->default_value("0.1"))
    ("ARATIO", "Radial scale length ratio for disk basis construction with doubleexpon",
     cxxopts::value<double>(ARATIO)->default_value("1.0"))
    ("HRATIO", "Vertical scale height ratio for disk basis construction with doubleexpon",
     cxxopts::value<double>(HRATIO)->default_value("1.0"))
    ("NUMR", "Size of radial grid",
     cxxopts::value<int>(NUMR)->default_value("2000"))
    ("PPOW", "Power-law index for power-law disk profile",
     cxxopts::value<double>(PPower)->default_value("4.0"))
    ("DWEIGHT", "Ratio of second disk relative to the first disk for disk basis construction with double-exponential",
     cxxopts::value<double>(DWEIGHT)->default_value("1.0"))
    ("RTRUNC", "Maximum disk radius for erf truncation of EOF conditioning density",
     cxxopts::value<double>(RTRUNC)->default_value("0.1"))
    ("RWIDTH", "Width for erf truncation for EOF conditioning density (ignored if zero)",
     cxxopts::value<double>(RWIDTH)->default_value("0.0"))
    ("RFACTOR", "Disk radial scaling factor for spherical deprojection model",
     cxxopts::value<double>(RFACTOR)->default_value("1.0"))
    ("RNUM", "Number of radial knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(RNUM)->default_value("200"))
    ("PNUM", "Number of azimthal knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(PNUM)->default_value("1"))
    ("TNUM", "Number of cos(theta) knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(TNUM)->default_value("80"))
    ("CMAPR", "Radial coordinate mapping type for cylindrical grid  (0=none, 1=rational fct)",
     cxxopts::value<int>(CMAPR)->default_value("1"))
    ("CMAPZ", "Vertical coordinate mapping type for cylindrical grid  (0=none, 1=rational fct)",
     cxxopts::value<int>(CMAPZ)->default_value("1"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
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


  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << options.help() << std::endl << std::endl
		<< "Examples: " << std::endl
		<< "\t" << "Use parameters read from a config file in INI style"  << std::endl
		<< "\t" << av[0] << " --config=gendisk.config"  << std::endl << std::endl
		<< "\t" << "Generate a template YAML config file current defaults called <template.yaml>"  << std::endl
		<< "\t" << av[0] << " --template" << std::endl << std::endl
		<< "\t" << "Override a single parameter in a config file from the command line"  << std::endl
		<< "\t" << av[0] << "--LMAX=8 --config=my.config" << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the YAML config file
  //
  if (vm.count("config")) {
    try {
      vm = LoadConfig(options, config);
    } catch (cxxopts::OptionException& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }

  // Enable new YAML cache header
  //
  if (vm.count("newcache")) {
    EmpCylSL::NewCache = true;
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
  if (vm.count("mtype")) {
    if (mtype.compare("exponential")==0)
      EmpCylSL::mtype = EmpCylSL::Exponential;
    else if (mtype.compare("gaussian")==0)
      EmpCylSL::mtype = EmpCylSL::Gaussian;
    else if (mtype.compare("plummer")==0)
      EmpCylSL::mtype = EmpCylSL::Plummer;
    else if (mtype.compare("power")==0) {
      EmpCylSL::mtype = EmpCylSL::Power;
      EmpCylSL::PPOW  = PPower;
    } else {
      if (myid==0) std::cout << "No EmpCylSL EmpModel named <"
			     << mtype << ">, valid types are: "
			     << "Exponential, Gaussian, Plummer, Power "
			     << "(not case sensitive)" << std::endl;
      MPI_Finalize();
      return 0;
    }
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

    if (myid==0)		// Report DiskType
      std::cout << "DiskType is <" << dtype << ">" << std::endl;
  }
  catch (const std::out_of_range& err) {
    if (myid==0) {
      std::cout << "DiskType error in configuraton file" << std::endl;
      std::cout << "Valid options are: ";
      for (auto v : dtlookup) std::cout << v.first << " ";
      std::cout << std::endl;
    }
    MPI_Finalize();
    return 0;
  }


  //====================
  // OpenMP control
  //====================

#ifdef HAVE_OMP_H
  omp_set_num_threads(nthrds);
#pragma omp parallel
  if (myid==0) {
    int numthrd = omp_get_num_threads();
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

                                // Create expansion only if needed . . .
  std::shared_ptr<EmpCylSL> expandd;

  expandd = std::make_shared<EmpCylSL>(NMAX2, LMAX2, MMAX, NORDER, ASCALE, HSCALE, NODD, cachefile);

#ifdef DEBUG
   std::cout << "Process " << myid << ": "
	     << " rmin="   << EmpCylSL::RMIN
	     << " rmax="   << EmpCylSL::RMAX
	     << " a="      << ASCALE
	     << " h="      << HSCALE
	     << " nmax2="  << NMAX2
	     << " lmax2="  << LMAX2
	     << " mmax="   << MMAX
	     << " nordz="  << NORDER
	     << " noddz="  << NODD
	     << std::endl  << std::flush;
#endif


   // Use these user models to deproject for the EOF spherical basis
   //
   if (vm.count("deproject")) {
     // The scale in EmpCylSL is assumed to be 1 so we compute the
     // height relative to the length
     //
     double H = HSCALE/ASCALE;

     // The model instance (you can add others in DiskModels.H).
     // It's MN or Exponential if not MN.
     //
     EmpCylSL::AxiDiskPtr model;

     if (dmodel.compare("MN")==0) // Miyamoto-Nagai
       model = std::make_shared<MNdisk>(1.0, H);
     else			// Default to exponential
       model = std::make_shared<Exponential>(1.0, H);

     if (RWIDTH>0.0) {
       model = std::make_shared<Truncated>(RTRUNC/ASCALE,
					   RWIDTH/ASCALE,
					   model);
       if (myid==0)
	 std::cout << "Made truncated model with R=" << RTRUNC/ASCALE
		   << " and W=" << RWIDTH/ASCALE << std::endl;
     }

     expandd->create_deprojection(H, RFACTOR, NUMR, RNUM, model);
   }

   // Regenerate EOF from analytic density
   //
   expandd->generate_eof(RNUM, PNUM, TNUM, dcond);

   // Basis orthgonality check
   //
   if (vm.count("ortho")) {
     std::ofstream out(runtag + ".ortho_check");
     expandd->ortho_check(out);
   }

  //===========================================================================
  // Shutdown MPI
  //===========================================================================

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
