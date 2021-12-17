/*
  Generates a Monte Carlo realization of a halo with an embedded
  disk using Jeans' equations.

  Assumptions:

  1) Spherical halo supplied in an input table

  2) Axisymmetric (but 3-dimensional) exponential disk with sech^2(Z/z)
     vertical profile

  3) Halo as spherical velocity ellipsoid

  4) Disk as axisymmetric velocity ellipsoid in the plane (that is,
     $\sigma_r = \sigma_\phi$ and $\sigma_z$ determined by solving
     Jeans' equations in cylindrical coordinates.

 Loosely based on Enrico Vesperini's initial.cc and diskANDhalo.cc
 (mostly rewritten)

 Added the basis expansion of the disk: 12/10/01. KHB

 Rewritten and debugged by MDW between 12/28/01-12/31/01.  

        Added command line parsing.  "gendisk -h" will list parameters.

        I removed everything but the "disk and halo" case.  Removed
        multiple inheritance.  Changed interface; expansions registered
        with DiskHalo class on construction

        Switched from biortho classes to expansion classes from the
        EXP code.

        Uses a vector of particle structures rather than a Matrix to
        store an pass phase space.

        Rewrote the particle component foliation code using a more
        general algorithm.

        Solution to Jeans' equations are now computed in parallel and
        tabulated.  

        Particles are stored on local nodes and written to disk by
        master.  

        Removed lots of other cruft.

 More debugging 03/05 by MDW

        Repaired EmpCylSL scaling

	Added additional debugging output

	Compared against expsl routines

	Removed orphaned parameters

 Updated to include gas disk using local Euler solution 04/08 by MDW

 Both constant scale height and isothermal gas disks 08/08 by MDW

 Multimass gas disk 11/08 by MDW

 Updates for constructing disk velocities from an evolved halo 01/20 by MDW

 Added double-exponential disk preconditioning 08/21 by MDW

*/
                                // C++/STL headers
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

#include <config.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

// EXP classes
//
#include <global.H>
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
#include <SphericalSL.H>
#include <DiskHalo.H>
#include <localmpi.H>


// Hydrogen fraction
//
const double f_H = 0.76;


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

  int          LMAX;
  int          NMAX;
  int          NUMR;
  int          SCMAP;
  int          NUMDF;
  double       RMIN;
  double       RCYLMIN;
  double       RCYLMAX;
  double       SCSPH;
  double       RSPHSL;
  double       DMFAC;
  double       RFACTOR;
  double       X0;
  double       Y0;
  double       Z0;
  double       U0;
  double       V0;
  double       W0;
  int          RNUM;
  int          PNUM;
  int          TNUM;
  int          VFLAG;
  int          DFLAG;
  bool         expcond;
  bool         LOGR;
  bool         CHEBY;
  int          CMAPR;
  int          CMAPZ;
  int          NCHEB;
  int          TCHEB;
  int          CMTYPE;
  int          NDR;
  int          NDZ;
  int          NHR;
  int          NHT;
  int          NDP;
  double       SHFAC;
  int          NMAX2;
  int          LMAX2;
  int          MMAX;
  int          NUMX;
  int          NUMY;
  int          NOUT;
  int          NORDER;
  int          NORDER1;
  int          NODD;
  bool         SELECT;
  bool         DUMPCOEF;
  int          DIVERGE;
  double       DIVERGE_RFAC;
  int          DIVERGE2;
  double       DIVERGE_RFAC2;
  int          DF;
  double       PPower;
  double       R_DF;
  double       DR_DF;
  double       Hratio;
  double       scale_height;
  double       scale_length;
  double       scale_lenfkN;
  double       disk_mass;
  double       gas_mass;
  double       gscal_length;
  double       ToomreQ;
  double       Temp;
  double       Tmin;
  //double       gen_ecut; // unused?
  bool         const_height=true;
  bool         images=false;
  bool         multi=false;
  bool         SVD;
  int          SEED;
  int          itmax;
  bool         DENS;
  bool         basis;
  bool         zero;
  bool         report;
  bool         ignore;
  bool         evolved;
  int          nhalo;
  int          ndisk;
  int          ngas;
  int          ngparam;
  string       hbods;
  string       dbods;
  string       gbods;
  string       suffix;
  string       centerfile="";
  string       halofile1;
  string       halofile2;
  string       cachefile;
  string       config;
  string       gentype;
  string       dtype;
  string       dmodel;
  string       mtype;
  string       ctype;
  
  const std::string mesg("Generates a Monte Carlo realization of a halo with an\n embedded disk using Jeans' equations\n");

  cxxopts::Options options(av[0], mesg);

  options.add_options()
    ("h,help", "Print this help message")
    ("T,template", "Write template options file with current and all default values",
     cxxopts::value<string>(config))
    ("f,input", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("deproject", "The EmpCylSL deprojection from specified disk model (EXP or MN)",
     cxxopts::value<string>(dmodel)->default_value("EXP"))
    ("hbods", "The output bodyfile for the halo",
     cxxopts::value<string>(hbods)->default_value("halo.bods"))
    ("gbods", "The output bodyfile for the gas disc",
     cxxopts::value<string>(gbods)->default_value("gas.bods"))
    ("dbods", "The output bodyfile for the stellar disc",
     cxxopts::value<string>(dbods)->default_value("disk.bods"))
    ("cachefile", "The cache file for the cylindrical basis",
     cxxopts::value<string>(cachefile)->default_value(".eof.cache.file"))
    ("ctype", "DiskHalo radial coordinate scaling type (one of: Linear, Log,Rat)",
     cxxopts::value<string>(ctype)->default_value("Log"))
    ("LMAX", "",
     cxxopts::value<int>(LMAX)->default_value("6"))
    ("NMAX", "",
     cxxopts::value<int>(NMAX)->default_value("20"))
    ("LMAX2", "",
     cxxopts::value<int>(LMAX2)->default_value("48"))
    ("NMAX2", "",
     cxxopts::value<int>(NMAX2)->default_value("48"))
    ("MMAX", "",
     cxxopts::value<int>(MMAX)->default_value("6"))
    ("NUMX", "",
     cxxopts::value<int>(NUMX)->default_value("256"))
    ("NUMY", "",
     cxxopts::value<int>(NUMY)->default_value("128"))
    ("DIVERGE", "",
     cxxopts::value<int>(DIVERGE)->default_value("0"))
    ("DIVERGE_RFAC", "",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("1.0"))
    ("DIVERGE2", "",
     cxxopts::value<int>(DIVERGE2)->default_value("0"))
    ("DIVERGE_RFAC2", "",
     cxxopts::value<double>(DIVERGE_RFAC2)->default_value("1.0"))
    ("nhalo", "Number of halo particles",
     cxxopts::value<int>(nhalo)->default_value("100000"))
    ("ndisk", "Number of disk particles",
     cxxopts::value<int>(ndisk)->default_value("100000"))
    ("ngas", "Number of gas disc particles",
     cxxopts::value<int>(ngas)->default_value("0"))
    ("ngparam", "Number of gas parameters",
     cxxopts::value<int>(ngparam)->default_value("0"))
    ("dtype", "",
     cxxopts::value<string>(dtype)->default_value("exponential"))
    ("NOUT", "",
     cxxopts::value<int>(NOUT)->default_value("18"))
    ("NODD", "",
     cxxopts::value<int>(NODD)->default_value("6"))
    ("NORDER", "",
     cxxopts::value<int>(NORDER)->default_value("1000"))
    ("NORDER1", "",
     cxxopts::value<int>(NORDER1)->default_value("10000"))
    ("NUMDF", "",
     cxxopts::value<int>(NUMDF)->default_value("1000"))
    ("VFLAG", "",
     cxxopts::value<int>(VFLAG)->default_value("31"))
    ("DFLAG", "",
     cxxopts::value<int>(DFLAG)->default_value("31"))
    ("expcond", "",
     cxxopts::value<bool>(expcond)->default_value("true"))
    ("report", "",
     cxxopts::value<bool>(report)->default_value("false"))
    ("ignore", "",
     cxxopts::value<bool>(ignore)->default_value("true"))
    ("evolved", "",
     cxxopts::value<bool>(evolved)->default_value("false"))
    ("multi", "Turn on multimass halo generation",
     cxxopts::value<bool>(multi)->default_value("false"))
    ("halofile1", "Halo spherical model table",
     cxxopts::value<string>(halofile1)->default_value("SLGridSph.model"))
    ("halofile2", "Number density profile for multimass halo generation",
     cxxopts::value<string>(halofile2)->default_value("SLGridSph.fake"))
    ("gentype", "",
     cxxopts::value<string>(gentype)->default_value("Asymmetric"))
    ("LOGR", "",
     cxxopts::value<bool>(LOGR)->default_value("true"))
    ("CHEBY", "",
     cxxopts::value<bool>(CHEBY)->default_value("false"))
    ("DENS", "",
     cxxopts::value<bool>(DENS)->default_value("true"))
    ("basis", "",
     cxxopts::value<bool>(basis)->default_value("false"))
    ("zero", "",
     cxxopts::value<bool>(zero)->default_value("true"))
    ("NCHEB", "",
     cxxopts::value<int>(NCHEB)->default_value("12"))
    ("TCHEB", "",
     cxxopts::value<int>(TCHEB)->default_value("12"))
    ("NDR", "",
     cxxopts::value<int>(NDR)->default_value("800"))
    ("SEED", "Random number seed for realization",
     cxxopts::value<int>(SEED)->default_value("11"))
    ("NDZ", "",
     cxxopts::value<int>(NDZ)->default_value("800"))
    ("NHR", "",
     cxxopts::value<int>(NHR)->default_value("800"))
    ("NHT", "",
     cxxopts::value<int>(NHT)->default_value("800"))
    ("NDP", "",
     cxxopts::value<int>(NDP)->default_value("800"))
    ("NUMR", "Size of radial grid for Spherical SL",
     cxxopts::value<int>(NUMR)->default_value("2000"))
    ("SHFAC", "",
     cxxopts::value<double>(SHFAC)->default_value("16"))
    ("disk_mass", "",
     cxxopts::value<double>(disk_mass)->default_value("0.05"))
    ("gas_mass", "",
     cxxopts::value<double>(gas_mass)->default_value("0.0"))
    ("scale_length", "",
     cxxopts::value<double>(scale_length)->default_value("0.01"))
    ("scale_height", "",
     cxxopts::value<double>(scale_height)->default_value("0.001"))
    ("ToomreQ", "",
     cxxopts::value<double>(ToomreQ)->default_value("1."))
    ("RMIN", "Minimum halo radius",
     cxxopts::value<double>(RMIN)->default_value("0.005"))
    ("RCYLMIN", "Minimum disk radius",
     cxxopts::value<double>(RCYLMIN)->default_value("0.001"))
    ("RCYLMAX", "Maximum disk radius",
     cxxopts::value<double>(RCYLMAX)->default_value("20.0"))
    ("SCMAP", "Turn on Spherical SL coordinate mapping (1, 2, 0=off)",
     cxxopts::value<int>(SCMAP)->default_value("1"))
    ("SCSPH", "Scale for Spherical SL coordinate mapping",
     cxxopts::value<double>(SCSPH)->default_value("1.0"))
    ("RSPHSL", "Maximum halo expansion radius",
     cxxopts::value<double>(RSPHSL)->default_value("47.5"))
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
    ("DWEIGHT", "Ratio of second disk relative to the first disk for disk basis construction with double-exponential",
     cxxopts::value<double>(DWEIGHT)->default_value("1.0"))
    ("RTRUNC", "Maximum disk radius for erf truncation of EOF conditioning density",
     cxxopts::value<double>(RTRUNC)->default_value("0.1"))
    ("RWIDTH", "Width for erf truncationofr EOF conditioning density (ignored if zero)",
     cxxopts::value<double>(RWIDTH)->default_value("0.0"))
    ("DMFAC", "Disk mass scaling factor for spherical deprojection model",
     cxxopts::value<double>(DMFAC)->default_value("1.0"))
    ("RFACTOR", "Disk radial scaling factor for spherical deprojection model",
     cxxopts::value<double>(RFACTOR)->default_value("1.0"))
    ("X0", "Disk-Halo x center position",
     cxxopts::value<double>(X0)->default_value("0.0"))
    ("Y0", "Disk-Halo y center position",
     cxxopts::value<double>(Y0)->default_value("0.0"))
    ("Z0", "Disk-Halo z center position",
     cxxopts::value<double>(Z0)->default_value("0.0"))
    ("U0", "Disk-Halo x velocity center position",
     cxxopts::value<double>(U0)->default_value("0.0"))
    ("V0", "Disk-Halo y velocity center position",
     cxxopts::value<double>(V0)->default_value("0.0"))
    ("W0", "Disk-Halo z velocity center position",
     cxxopts::value<double>(W0)->default_value("0.0"))
    ("RNUM", "Number of radial knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(RNUM)->default_value("200"))
    ("PNUM", "Number of azimthal knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(PNUM)->default_value("1"))
    ("TNUM", "Number of cos(theta) knots for EmpCylSL basis construction quadrature",
     cxxopts::value<int>(TNUM)->default_value("80"))
    ("CMAPR", "Radial coordinate mapping type for cylindrical grid  (0=none, 1=rational fct)",
     cxxopts::value<int>(CMAPR)->default_value("1"))
    ("CMAPZ", "Verticall coordinate mapping type for cylindrical grid  (0=none, 1=rational fct)",
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
    // Do not overwrite existing config file
    //
    if (std::filesystem::exists(config)) {
      if (myid == 0)
	std::cerr << av[0] << ": config file <" << config
		  << "> exists, will not overwrite" << std::endl;
      MPI_Finalize();
      return 0;
    }

    NOUT = std::min<int>(NOUT, NORDER);

    // Write template file
    //
    if (myid==0) SaveConfig(vm, options, config);

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
		<< "\t" << av[0] << " --input=gendisk.config"  << std::endl << std::endl
		<< "\t" << "Generate a template config file in INI style from current defaults"  << std::endl
		<< "\t" << av[0] << " --conf=template.config" << std::endl << std::endl
		<< "\t" << "Override a single parameter in a config file from the command line"  << std::endl
		<< "\t" << av[0] << "--LMAX=8 --conf=template.config" << std::endl << std::endl;
    }
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
  
  // Enable new YAML cache header
  //
  if (vm.count("newcache")) {
    EmpCylSL::NewCache = true;
  }

  if (vm.count("spline")) {
    SphericalModelTable::linear = 0;
  } else {
    SphericalModelTable::linear = 1;
  }

  // Set EmpCylSL mtype.  This is the spherical function used to
  // generate the EOF basis.  If "deproject" is set, this will be
  // overriden in EmpCylSL.
  //

				// Convert mtype string to lower case
  std::transform(mtype.begin(), mtype.end(), mtype.begin(),
		 [](unsigned char c){ return std::tolower(c); });

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

  // Set DiskType.  This is the functional form for the disk used to
  // condition the basis.
  //
				// Convert dtype string to lower case
  std::transform(dtype.begin(), dtype.end(), dtype.begin(),
		 [](unsigned char c){ return std::tolower(c); });

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
  // Cheb order for
  // SphericalModel DF 
  //====================

  SphericalModelTable::chebyN = TCHEB;

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
  // Random variates
  //====================

  random_gen.seed(SEED);
  std::uniform_real_distribution<> unit;
  std::normal_distribution<> unitN;

  //====================
  // Okay, now begin ...
  //====================

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  // set_fpu_handler();         // Make gdb trap FPU exceptions
  set_fpu_gdb_handler();	// Make gdb trap FPU exceptions
#endif
  
  int n_particlesH, n_particlesD, n_particlesG;
  
  if (suffix.size()>0) {
    hbods = hbods + "." + suffix;
    dbods = dbods + "." + suffix;
    gbods = gbods + "." + suffix;
  }
  
  // Divvy up the particles by core
  //
  n_particlesH = nhalo/numprocs;
  if (myid==0) n_particlesH = nhalo - n_particlesH*(numprocs-1);
  
  n_particlesD = ndisk/numprocs;
  if (myid==0) n_particlesD = ndisk - n_particlesD*(numprocs-1);
  
  n_particlesG = ngas/numprocs;
  if (myid==0) n_particlesG = ngas  - n_particlesG*(numprocs-1);
  
  
#ifdef DEBUG  
  std::cout << "Processor " << myid << ": n_particlesH=" << n_particlesH
	    << std::endl
	    << "Processor " << myid << ": n_particlesD=" << n_particlesD
	    << std::endl
	    << "Processor " << myid << ": n_particlesG=" << n_particlesG
	    << std::endl;
#endif
  
  if (nhalo + ndisk + ngas <= 0) {
    if (myid==0) std::cout << "You have specified zero particles!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
    exit(0);
  }
  
  // Vectors to contain phase space Particle structure is defined in
  // Particle.H
  //
  vector<Particle> dparticles, hparticles;
  
  //
  // Disk halo grid parameters
  //
  //      Prevent model evaulation inside of either grid---------+
  //                                       |                     |
  //                                       V                     V
  DiskHalo::RDMIN       = std::max<double>(RCYLMIN*scale_length, RMIN);
  DiskHalo::RHMIN       = RMIN;
  DiskHalo::RHMAX       = RSPHSL;
  DiskHalo::RDMAX       = RCYLMAX*scale_length;
  DiskHalo::NDR         = NDR;
  DiskHalo::NDZ         = NDZ;
  DiskHalo::NHR         = NHR;
  DiskHalo::NHT         = NHT;
  DiskHalo::NDP         = NDP;
  DiskHalo::SHFACTOR    = SHFAC;
  DiskHalo::COMPRESSION = DMFAC;
  DiskHalo::NUMDF       = NUMDF;
  DiskHalo::Q           = ToomreQ;
  DiskHalo::R_DF        = R_DF;
  DiskHalo::DR_DF       = DR_DF;
  DiskHalo::SEED        = SEED;
  DiskHalo::VFLAG       = static_cast<unsigned int>(DFLAG);
  DiskHalo::CHEBY       = CHEBY;
  DiskHalo::NCHEB       = NCHEB;

  if (vm.count("itmax"))  DiskHalo::ITMAX    = itmax;
  if (vm.count("allow"))  DiskHalo::ALLOW    = true;
  if (vm.count("nomono")) DiskHalo::use_mono = false;
  if (suffix.size())      DiskHalo::RUNTAG   = suffix;
  
  AddDisk::use_mpi      = true;
  AddDisk::Rmin         = RMIN;
  
  //===========================Spherical expansion=============================
  
  // SLGridSph::diverge = DIVERGE;
  // SLGridSph::divergexp = DIVERGE_RFAC;
  SLGridSph::model_file_name = halofile1;
  
  SphericalSL::RMIN = RMIN;
  SphericalSL::RMAX = RSPHSL;
  SphericalSL::NUMR = NUMR;
  // Create expansion only if needed . . .
  std::shared_ptr<SphericalSL> expandh;
  if (nhalo) {
    expandh = std::make_shared<SphericalSL>(nthrds, LMAX, NMAX, SCMAP, SCSPH);
  }

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
  EmpCylSL::USESVD      = SVD;
  EmpCylSL::PCAVAR      = SELECT;
  EmpCylSL::CACHEFILE   = cachefile;

  if (basis) EmpCylSL::DENS = true;

                                // Create expansion only if needed . . .
  std::shared_ptr<EmpCylSL> expandd;
  bool save_eof = false;

  if (ndisk) {

    expandd = std::make_shared<EmpCylSL>(NMAX2, LMAX2, MMAX, NORDER, ASCALE, HSCALE, NODD);

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
	return 0;
      }
    }

    // Use these user models to deproject for the EOF spherical basis
    //
    if (vm.count("deproject")) {
      // The scale in EmpCylSL is assumed to be 1 so we compute the
      // height relative to the length
      //
      double H = scale_height/scale_length;

      // The model instance (you can add others in DiskModels.H).
      // It's MN or Exponential if not MN.
      //
      EmpCylSL::AxiDiskPtr model;

      if (dmodel.compare("MN")==0) // Miyamoto-Nagai
	model = std::make_shared<MNdisk>(1.0, H);
      else			// Default to exponential
	model = std::make_shared<Exponential>(1.0, H);

      if (RWIDTH>0.0) {
	model = std::make_shared<Truncated>(RTRUNC/scale_length,
					      RWIDTH/scale_length,
					      model);
	if (myid==0)
	  std::cout << "Made truncated model with R=" << RTRUNC/scale_length
		    << " and W=" << RWIDTH/scale_length << std::endl;
      }

      expandd->create_deprojection(H, RFACTOR, NUMR, RNUM, model);
    }

    // Regenerate EOF from analytic density
    //
    if (expcond and not save_eof) {
      expandd->generate_eof(RNUM, PNUM, TNUM, dcond);
      save_eof = true;
    }

    // Basis orthgonality check
    //
    if (vm.count("ortho")) {
      std::ofstream out(runtag + ".ortho_check");
      expandd->ortho_check(out);
    }
  }

  //====================Create the disk & halo model===========================

  std::shared_ptr<DiskHalo> diskhalo;

  if (multi) {
    if (myid==0) std::cout << "Initializing for a MULTI-MASS halo . . . " << std::flush;
    diskhalo =
      std::make_shared<DiskHalo>
      (expandh, expandd,
       scale_height, scale_length, disk_mass, ctype,
       halofile1, DIVERGE,  DIVERGE_RFAC,
       halofile2, DIVERGE2, DIVERGE_RFAC2,
       DiskHalo::getDiskGenType[gentype]);
    if (myid==0) std::cout << "done" << std::endl;

  } else {

    if (myid==0) std::cout << "Initializing for a SINGLE-MASS halo . . . " << std::flush;
    diskhalo = std::make_shared<DiskHalo>
      (expandh, expandd,
       scale_height, scale_length, 
       disk_mass, ctype, halofile1,
       DF, DIVERGE, DIVERGE_RFAC,
       DiskHalo::getDiskGenType[gentype]);
    if (myid==0) std::cout << "done" << std::endl;
  }
  
  std::ifstream center(centerfile);
  if (center) {

    bool ok = true;

    center >> X0;
    if (center.fail()) ok = false;

    center >> Y0;
    if (center.fail()) ok = false;

    center >> Z0;
    if (center.fail()) ok = false;

    if (ok) {
      diskhalo->set_pos_origin(X0, Y0, Z0);
      if (myid==0) std::cout << "Using position origin: " 
			     << X0 << ", " << Y0 << ", " << Z0 << std::endl;
    }

    center >> U0;
    if (center.fail()) ok = false;

    center >> V0;
    if (center.fail()) ok = false;

    center >> W0;
    if (center.fail()) ok = false;

    if (ok) {
      diskhalo->set_vel_origin(U0, V0, W0);
      if (myid==0) std::cout << "Using velocity origin: " 
			     << U0 << ", " << V0 << ", " << W0 << std::endl;
    }
  }

                                // Make zero center of mass and
                                // center of velocity
  diskhalo->zero_com(zero);
  diskhalo->zero_cov(zero);
  
  //===========================================================================

                                // Open output file (make sure it exists
                                // before realizing a large phase space)
  std::ofstream out_halo, out_disk;
  if (myid==0) {
    if (not evolved and n_particlesH) {
      out_halo.open(hbods);
      if (!out_halo) {
	cout << "Could not open <" << hbods << "> for output\n";
	MPI_Abort(MPI_COMM_WORLD, 4);
	exit(0);
      }
    }

    if (ndisk) {
      out_disk.open(dbods);
      if (!out_disk) {
	std::cout << "Could not open <" << dbods << "> for output" << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 4);
	exit(0);
      }
    }
  }

  //=================Make the phase space coordinates==========================

  if (evolved) {		// ---------------------------
				// Use existing halo body file
    std::ifstream hin(hbods);	// ---------------------------
    
    if (hin) {
      std::string line;
      std::getline(hin, line);
      std::istringstream sin(line);

      int niatr, ndatr;
      sin >> nhalo;
      sin >> niatr;
      sin >> ndatr;
      
      // Divvy up the particles by core.  The root node gets any
      // remainder.
      //
      n_particlesH = nhalo/numprocs;

      int ibeg = 0;
      int iend = nhalo - n_particlesH*(numprocs-myid-1);
      
      if (myid>0) {
	ibeg = nhalo - n_particlesH*(numprocs-myid);
	for (int i=0; i<ibeg; i++) std::getline(hin, line);
      }

      Particle P(niatr, ndatr);

      for (int i=ibeg; i<iend; i++) {
	std::getline(hin, line);
	std::istringstream sin(line);
	sin >> P.mass;
	for (int k=0; k<3; k++)     sin >> P.pos[k];
	for (int k=0; k<3; k++)     sin >> P.vel[k];
	for (int k=0; k<niatr; k++) sin >> P.iattrib[k];
	for (int k=0; k<ndatr; k++) sin >> P.dattrib[k];
	hparticles.push_back(P);
      }

    } else {
      // Error message
      if (myid==0)
	std::cout << "Could not read halo file <" << hbods
		  << "> . . . quitting" << std::endl;

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(-1);
    }

    std::cout << "Process " << myid << " has " << hparticles.size()
	      << " halo particles" << std::endl;

    if (myid==0)
      std::cout << "Generating halo tables for input halo . . . "
		<< std::flush;

    if (multi) {
      diskhalo->set_halo_table_multi(hparticles);
    } else {
      diskhalo->set_halo_table_single(hparticles);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) std::cout << "done" << std::endl;

  } else {			// ---------------------------
				// Generate new halo body file
    if (nhalo) {		// ---------------------------
      if (multi) {
	if (myid==0) std::cout << "Generating halo phase space . . . " << std::flush;
	diskhalo->set_halo(hparticles, nhalo, n_particlesH);
      } else {
	if (myid==0) std::cout << "Generating halo coordinates . . . " << std::flush;
	diskhalo->set_halo_coordinates(hparticles, nhalo, n_particlesH);
	MPI_Barrier(MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==0) std::cout << "done" << std::endl;
    }
  }

  if (ndisk) {
    if (myid==0) std::cout << "Generating disk coordinates . . . " << std::flush;
    diskhalo->set_disk_coordinates(dparticles, ndisk, n_particlesD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) std::cout << "done" << std::endl;
  }

  if (nhalo) {
    if (myid==0) std::cout << "Beginning halo accumulation . . . " << std::flush;
    expandh->accumulate(hparticles);
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << "done" << std::endl;
      if (DUMPCOEF) {
	std::cout << "Dumping coefficients halo . . . " << std::flush;
	ostringstream sout;
	sout << "halo_coefs.";
	if (suffix.size()>0) sout << suffix;
	else                 sout << "dump";
	ofstream out(sout.str());
	if (out) expandh->dump_coefs(out, false);
	std::cout << "done" << std::endl;
      }
      if (vm.count("probe")) {
	std::cout << "Dumping a probe through the halo . . . " << std::flush;
	ostringstream sout;
	sout << "halo_probe.";
	if (suffix.size()>0) sout << suffix;
	else                 sout << "dump";
	ofstream out(sout.str());
	if (out) {
	  // Header
	  out << "#" << std::setw(15) << std::right << "radius |"
	      << std::setw(16) << std::right << "density |"
	      << std::setw(16) << std::right << "potential |"
	      << std::setw(16) << std::right << "d(pot)/dr |"
	      << std::setw(16) << std::right << "d(pot)/dt |"
	      << std::setw(16) << std::right << "d(pot)/dp |"
	      << std::endl
	      << "#" << std::setw(15) << std::right << "[1] |"
	      << std::setw(16) << std::right << "[2] |"
	      << std::setw(16) << std::right << "[3] |"
	      << std::setw(16) << std::right << "[4] |"
	      << std::setw(16) << std::right << "[5] |"
	      << std::setw(16) << std::right << "[6] |"
	      << std::endl;

	  // Do the probe
	  //
	  const int nn = 1000;
	  double dr = (log(RSPHSL) - log(RMIN))/(nn-1);
	  for (int n=0; n<nn; n++) {
	    double rr = RMIN*exp(dr*n);
	    double dens, potl, dpr, dpt, dpp;
	    expandh->determine_fields_at_point(rr, 0.0, 0.0,
					       &dens, &potl, &dpr, &dpt, &dpp);
	    out << std::setw(16) << rr
		<< std::setw(16) << dens
		<< std::setw(16) << potl
		<< std::setw(16) << dpr
		<< std::setw(16) << dpt
		<< std::setw(16) << dpp
		<< std::endl;
	  }
	} else {
	  std::cout << "Could not open file: " << sout.str() << std::endl;
	}
	std::cout << "done" << std::endl;
      }
    }

  }
  
  if (ndisk) {
    if (myid==0) std::cout << "Beginning disk accumulation . . . " << std::flush;
    expandd->setup_accumulation();

    if (!save_eof and !expcond) {
      expandd->setup_eof();
      if (nthrds>1)
	expandd->accumulate_eof_thread(dparticles, report);
      else
	expandd->accumulate_eof(dparticles, report);
      MPI_Barrier(MPI_COMM_WORLD);

      if (myid==0) std::cout << "done" << std::endl;
  
      if (myid==0) std::cout << "Making the EOF . . . " << std::flush;
      expandd->make_eof();
      MPI_Barrier(MPI_COMM_WORLD);
    }

    if (myid==0) {
      std::cout << "done" << std::endl;
      std::cout << "Making disk coefficients . . . " << std::flush;
    }

    expandd->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << "done" << std::endl;
      std::cout << "Reexpand . . . " << std::flush;
    }

    if (nthrds>1)
      expandd->accumulate_thread(dparticles, 0, report);
    else
      expandd->accumulate(dparticles, 0, report);

    expandd->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      std::cout << "done" << std::endl;
      if (DUMPCOEF) {
	std::cout << "Dumping coefficients . . . " << std::flush;
	ostringstream sout;
	sout << "disk_coefs.";
	if (suffix.size()>0) sout << suffix;
	else                 sout << "dump";
	ofstream out(sout.str());
	if (out) expandd->dump_coefs(out);
	std::cout << "done" << std::endl;
      }
    }

    if (NORDER1<NORDER) {
      if (myid==0) std::cout << "Restricting order from " << NORDER 
			     << " to " << NORDER1 << " . . . " << std::flush;
      expandd->restrict_order(NORDER1);
      if (myid==0) std::cout << "done" << std::endl;
    }

    if (images && myid==0) {
      std::cout << "Images . . . " << std::flush;
      std::ostringstream dumpname;
      dumpname << "images.0";
      expandd->dump_images(dumpname.str(), 5.0*ASCALE, 5.0*HSCALE, 64, 64, true);
      std::cout << "done" << std::endl;
    }
  }
  

  //===========================Diagnostics=====================================

                                // For examining the coverage, etc.
                                // Images can be contoured in SM using
  if (myid==0 && basis) {	// the "ch" file type
    
    std::cout << "Dumping basis images . . . " << std::flush;
    
    if (ndisk) {
      int nout = 200;
      string dumpstr = runtag + ".dump";
      expandd->dump_basis(dumpstr, 0);
      expandd->dump_images(runtag, 5.0*scale_length, 5.0*scale_height,
			   nout, nout, false);
      expandd->dump_images_basis(runtag, 5.0*scale_length, 5.0*scale_height,
				 nout, nout, false, 0, MMAX, 0, NORDER-1);
    }


    if (nhalo) {
      string extn("test");
      expandh->dump_basis(extn);
    }
    
    if (nhalo) {
      
      const int nstr = 5;
      const char *names[nstr] = {".dens", ".potl", ".potr", ".pott", ".potp"};
      std::vector<std::ofstream> out(nstr);
      
      int nout = 200;
      double rmax = 6.0*scale_length;
      double x, y, dr = 2.0*rmax/(nout-1);
      float f;
    
      for (int i=0; i<nstr; i++) {
        string name("halo");
        name += names[i];
        out[i].open(name);
        
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
      }
      
      double r, theta, phi;
      double dens, potl, potr, pott, potp;
    
      for (int j=0; j<nout; j++) {
        y = -rmax + dr*j;
      
        for (int i=0; i<nout; i++) {
          x = -rmax + dr*i;
        
          r = sqrt(x*x + y*y);
          theta = 0.5*M_PI;
          phi = atan2(y, x);
        
          expandh->determine_fields_at_point(r, theta, phi,
                                             &dens, &potl, 
                                             &potr, &pott, &potp);
        
          out[0].write((char *)&(f=dens), sizeof(float));
          out[1].write((char *)&(f=potl), sizeof(float));
          out[2].write((char *)&(f=potr), sizeof(float));
          out[3].write((char *)&(f=pott), sizeof(float));
          out[4].write((char *)&(f=potp), sizeof(float));
        }
        
      }
    
      for (int i=0; i<nstr; i++) out[i].close();
    }

    if (ndisk) {

      const int nstr = 5;
      const char *names[nstr] = {".dens", ".pot", ".fr", ".fz", ".fp"};
      std::vector<std::ofstream> out(nstr);
    
      int nout = 200;
      double rmax = DiskHalo::RDMAX;
      double x, y, dr = 2.0*rmax/(nout-1);
      float f;
    
      for (int i=0; i<nstr; i++) {
        string name("disk");
        name += names[i];
        out[i].open(name);
        
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
      }
    
      double z = 0.0, d0, p0, d, p, fr, fz, fp;
    
      for (int j=0; j<nout; j++) {
        y = -rmax + dr*j;
      
        for (int i=0; i<nout; i++) {
          x = -rmax + dr*i;
        
	  if (x<0.0)
	    expandd->accumulated_eval(fabs(x), y, M_PI, p0, p, fr, fz, fp);
	  else
	    expandd->accumulated_eval(x, y,  0.0, p0, p, fr, fz, fp);


          d = expandd->accumulated_dens_eval(sqrt(x*x + y*y), z, atan2(y, x), d0);
        
          
          out[0].write((char *)&(f=d ), sizeof(float));
          out[1].write((char *)&(f=p ), sizeof(float));
          out[2].write((char *)&(f=fr), sizeof(float));
          out[3].write((char *)&(f=fz), sizeof(float));
          out[4].write((char *)&(f=fp), sizeof(float));
        }
        
      }
    
      for (int i=0; i<nstr; i++) out[i].close();
    }
    
    std::cout << "done" << std::endl;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  //====================Make the phase space velocities========================

  if (!multi and !evolved) {
    if (myid==0) std::cout << "Generating halo velocities . . . " << std::flush;
    diskhalo->set_vel_halo(hparticles);
    if (myid==0) std::cout << "done" << std::endl;
  }
  
  if (myid==0) std::cout << "Generating disk velocities . . . " << std::flush;
  diskhalo->set_vel_disk(dparticles);
  if (myid==0) std::cout << "done" << std::endl;
  

  //====================All done: write it out=================================

  if (not evolved) {
    if (myid==0) std::cout << "Writing phase space file for halo . . . " << std::flush;
    diskhalo->write_file(out_halo, hparticles);
    if (myid==0) std::cout << "done" << std::endl;
    out_halo.close();
  }

  if (myid==0) std::cout << "Writing phase space file for disk . . . " << std::flush;
  diskhalo->write_file(out_disk, dparticles);
  if (myid==0) std::cout << "done" << std::endl;
  out_disk.close();
                                // Diagnostic . . .
  diskhalo->virial_ratio(hparticles, dparticles);

  std::ofstream outprof("profile.diag");
  diskhalo->profile(outprof, dparticles, 3.0e-3*ASCALE, 5.0*ASCALE, 100);

  //====================Compute gas particles==================================

  if (myid==0 && n_particlesG) {
    std::cout << "Computing gas particles . . . " << std::endl;

    // UNITS
    // -------------------
				// cm
    const double pc = 3.08568025e18;
				// proton mass
    const double m_p = 1.67262158e-24;
				// g
    const double msun = 1.98892e33; //
				// cgs
    const double G = 6.67300e-08;
				// cgs
    const double boltz = 1.3806503e-16;

    double T = Temp;

    
    double Lunit = 3.0e5*pc;	// Virial radius
    double Munit = 1.0e12*msun;	// Virial mass
    double Tunit = sqrt(Lunit*Lunit*Lunit/(Munit*G));
    double Vunit = Lunit/Tunit;

    // Fac = kT*R_vir/(G*m_p*M_vir)
    // where M_vir = 1e12 Msun, R_vir=300 kpc
    //
    double fac = T/(G*m_p*Munit/(Lunit*boltz));

    // Thermal velocity in system units
    //
    double mm   = f_H*m_p + (1.0-f_H)*4.0*m_p;
    double vthermal = sqrt( (boltz*T)/mm ) / Vunit;
    double vmin2 = (boltz*Tmin/mm) / (Vunit*Vunit);

    // Adjust scale for multimass gas
    //
    double Scale_Length = gscal_length;
    if (scale_lenfkN > 0.0) gscal_length = scale_lenfkN;

    // Compute using Jeans theorem
    //
    double rmin = RMIN;
    double rmax = 10.0*gscal_length;
    double zmin = 0.001*scale_height;
    int   nrint = 200;
    int   nzint = 400;
    
    vector< vector<double> > zrho, zmas, vcir;
    double r, R, dR = (rmax - rmin)/(nrint-1);
    double z, dz = (log(rmax) - log(zmin))/(nzint-1);

    double p0, p, fr, fz, fp, dens, potl, potr, pott, potp;

    std::cout << "Const_height=" << (const_height ? "True" : "False") << std::endl;

    if (const_height) {

      for (int i=0; i<nrint; i++) {
	R = rmin + dR*i;

	vector<double> lrho(nzint), trho(nzint), tcir(nzint), tmas(nzint, 0);

	for (int j=0; j<nzint; j++) {
	  z = zmin*exp(dz*j);
	  r = sqrt(R*R + z*z);
	  
	  double pot=0.0, frt0=0.0, fzt0=0.0;
	  if (expandd) {
	    expandd->accumulated_eval(R, z, 0, p0, p, fr, fz, fp);
	    frt0 += -fr;
	    fzt0 += -fz;
	    pot += p0;
	  }
	  if (expandh) {
	    expandh->determine_fields_at_point(r, acos(z/(r+1.0e-8)), 0.0,
					       &dens, &potl, 
					       &potr, &pott, &potp);
	    
	    frt0 += potr;
	    fzt0 += (potr*z + pott*R*R/(r*r))/r;
	    pot += potl;
	  }
	  
	  trho[j] = fzt0*scale_height;
	  tcir[j] = sqrt(max<double>(R*frt0-R*trho[j]/Scale_Length, 0.0));
	}
	
	for (int j=0; j<nzint; j++) 
	  tmas[j] = 1.0 - exp(-zmin*exp(dz*j)/scale_height);
	
	zrho.push_back(trho);
	zmas.push_back(tmas);
	vcir.push_back(tcir);
      }

      //
      // Vertical table
      //
      string ztable("ztable.dat");
      std::cout << "Writing " << setw(15) << right << ztable
		<< " [gas] . . . " << std::flush;
      ofstream ztest(ztable);
      for (int i=0; i<nrint; i++) {
	for (int j=0; j<nzint; j++) {
	  ztest << setw(15) << rmin + dR*i
		<< setw(15) << zmin*exp(dz*j)
		<< setw(15) << zrho[i][j]
		<< setw(15) << zrho[i][j]*Vunit*Vunit*mm/boltz
		<< setw(15) << zmas[i][j]
		<< setw(15) << vcir[i][j]
		<< endl;
	}
	ztest << endl;
      }
      ztest.close();
      std::cout << "done" << std::endl;
      
    } else {

      for (int i=0; i<nrint; i++) {
	R = rmin + dR*i;



	vector<double> lrho(nzint), trho(nzint), tcir(nzint), tmas(nzint, 0);

	for (int j=0; j<nzint; j++) {
	  z = zmin*exp(dz*j);
	  r = sqrt(R*R + z*z);
	  
	  double frt0=0.0, fzt0=0.0;
	  if (expandd) {
	    expandd->accumulated_eval(R, z, 0, p0, p, fr, fz, fp);
	    frt0 += -fr;
	    fzt0 += -fz;
	  }
	  if (expandh) {
	    expandh->determine_fields_at_point(r, acos(z/(r+1.0e-8)), 0.0,
					       &dens, &potl, 
					       &potr, &pott, &potp);
	    frt0 += potr;
	    fzt0 += (potr*z + pott*R*R/(r*r))/r;
	  }
	  
	  trho[j] = -fzt0/(vthermal*vthermal);
	  tcir[j] = sqrt(max<double>(R*frt0-R*vthermal*vthermal/Scale_Length, 0.0));
	}
	
	double mass = 0.0;
	double zfac = 1.0 - exp(-dz);
				    
	lrho[0] = 0.0;
	for (int j=1; j<nzint; j++) 
	  lrho[j] = lrho[j-1] + 0.5*(trho[j-1] + trho[j]) * zmin*exp(dz*j)*zfac;
	
	for (int j=1; j<nzint; j++) 
	  tmas[j] = tmas[j-1] + 0.5*(exp(lrho[j-1]) + exp(lrho[j])) * zmin*exp(dz*j)*zfac;
	
	for (int j=0; j<nzint; j++) {
	  if (tmas[nzint-1]>0.0 && !std::isnan(tmas[nzint-1])) {
	    trho[j]  = exp(lrho[j])/tmas[nzint-1];
	    tmas[j] /= tmas[nzint-1];
	  } else {
	    trho[j] = 0.0;
	    if (j==0) tmas[j] = 0.0;
	    else      tmas[j] = 1.0;
	  }
	}
	zrho.push_back(trho);
	zmas.push_back(tmas);
	vcir.push_back(tcir);
      }


      //
      // Vertical table
      //
      std::cout << "Writing ztable.dat [gas] . . . " << std::flush;
      std::ofstream ztest("ztable.dat");
      for (int i=0; i<nrint; i++) {
	for (int j=0; j<nzint; j++) {
	  ztest << setw(15) << rmin + dR*i
		<< setw(15) << zmin*exp(dz*j)
		<< setw(15) << zrho[i][j]
		<< setw(15) << zmas[i][j]
		<< setw(15) << vcir[i][j]
		<< endl;
	}
	ztest << endl;
      }
      ztest.close();
      std::cout << "done" << std::endl;
      
    }

    // 
    // Prepare output stream
    //
    ofstream outps("gas.bods");
    if (!outps) {
      cerr << "Couldn't open <" << "gas.bods" << "> for output\n";
      exit (-1);
    }

    const int ITMAX=1000;
    const int NREPORT=1000;
    
    //
    // Maximum enclosed disk mass given rmax
    //
    double rmx2 = 1.5*rmax;
    double mmx2 = 1.0 - (1.0 + rmx2/gscal_length)*exp(-rmx2/gscal_length);
    double mmax = 1.0 - (1.0 + rmax/gscal_length)*exp(-rmax/gscal_length);
    double mfac = 1.0 - (1.0 + rmax/Scale_Length)*exp(-rmax/Scale_Length);

    //
    // Trimmed Gaussian
    //
    double minK=0.0, maxK=1.0, sigma = 3.0;
    if (sigma>0) {
      minK = 0.5*(1.0+erf(-0.5*sigma));
      maxK = 0.5*(1.0+erf( 0.5*sigma));
    }

    double gmass, gmass0 = gas_mass/ngas;
    double KE=0.0, VC=0.0;
    vector<double> mc2(nzint);

    gmass = gmass0;
    fr = fz = potr = 0.0;

    outps << setw(8) << ngas
	  << setw(6) << 0 << setw(6) << ngparam << endl;

    for (int n=0; n<ngas; n++) {

      double F, dF, M=mmax*unit(random_gen), Z=unit(random_gen);
      double R = M*rmax, phi=2.0*M_PI*unit(random_gen), x, y, z, rr, vc;
      double ax, ay, az;

				// Narrow with bisection
      double rm = 0.0, rp = rmx2;
      double fm = -M, fp = mmx2 - M;
      for (int j=0; j<15; j++) {
	R = 0.5*(rm + rp);
	F = 1.0 - M - (1.0 + R/gscal_length)*exp(-R/gscal_length);
	if (fm*F<0.0) {
	  rp = R;
	  fp = F;
	} else {
	  rm = R;
	  fm = F;
	}
      }
				// Polish with Newton-Raphson
      for (int j=0; j<ITMAX; j++) {
	F = 1.0 - M - (1.0 + R/gscal_length)*exp(-R/gscal_length);
	dF = R/(gscal_length*gscal_length)*exp(-R/gscal_length);
	R += -F/dF;
	if (fabs(F/dF)<1.0e-12) break;
      }
    
      int indr = static_cast<int>(floor(R/dR));
      if (indr<0) indr=0;
      if (indr>nrint-2) indr=nrint-2;
      double a = (dR*(indr+1) - R)/dR;
      double b = (R - indr*dR)/dR;

      vector<double> mz(nzint), vz(nzint);
      for (int j=0; j<nzint; j++) {
	mz[j] = a*zmas[indr][j] + b*zmas[indr+1][j];
	vz[j] = a*vcir[indr][j] + b*vcir[indr+1][j];
      }
      for (int j=0; j<nzint; j++) mz[j] /= mz[nzint-1];
      
      if (const_height) {
	for (int j=0; j<nzint; j++) 
	  mc2[j] = max<double>(a*zrho[indr][j] + b*zrho[indr+1][j], vmin2);
      }

      int indz = max<int>(0, min<int>(nzint-2, Vlocate(Z, mz)));

      a = (mz[indz+1] - Z)/(mz[indz+1] - mz[indz]);
      b = (Z - mz[indz  ])/(mz[indz+1] - mz[indz]);

      vc = fabs(a*vz[indr] + b*vz[indr+1]);

      z = zmin*exp(dz*(a*indz + b*(indz+1)));
      if (unit(random_gen)<0.5) z *= -1.0;
      rr = sqrt(R*R + z*z);

      if (const_height) {
	vthermal = a*mc2[indz] + b*mc2[indz+1];
	vthermal = sqrt(max<double>(vmin2, vthermal));
      }

      double sinp = sin(phi), cosp = cos(phi);
      x = R*cosp;
      y = R*sinp;

      double u = -vc*sinp + vthermal*norminv(unitN(random_gen));
      double v =  vc*cosp + vthermal*norminv(unitN(random_gen));
      double w =  vthermal*norminv(unitN(random_gen));
      
      gmass = gmass0*exp(-R*(1.0/Scale_Length - 1.0/gscal_length)) * 
	mmax*gscal_length*gscal_length/(mfac*Scale_Length*Scale_Length);

      outps << setw(18) << gmass
	    << setw(18) << R*cos(phi)
	    << setw(18) << R*sin(phi)
	    << setw(18) << z
	    << setw(18) << u
	    << setw(18) << v
	    << setw(18) << w;
      for (int k=0; k<ngparam; k++) outps << setw(18) << 0.0;
      outps << endl;
    
      if (expandd)
	expandd->accumulated_eval(R, z, phi, p0, p, fr, fz, fp);

      if (expandh)
	expandh->determine_fields_at_point(rr, acos(z/(rr+1.0e-8)), 0.0,
					   &dens, &potl, 
					   &potr, &pott, &potp);
      KE += 0.5*gmass*(u*u + v*v + w*w);

      VC += gmass*(-rr*potr + R*fr + z*fz);

      if (!((n+1)%NREPORT)) std::cout << "\r." << n+1 << std::flush;
    }

    std::cout << endl << "Done!" << std::endl;

    std::cout << "****************************" << std::endl
	      << "  Gas disk"                   << std::endl
	      << "----------------------------" << std::endl
	      << "  KE       = " << KE << std::endl
	      << "  VC       = " << VC << std::endl;
    if (VC<0.0)
      std::cout << " -2T/W     = " << -2.0*KE/VC << std::endl;
    std::cout << "****************************"  << std::endl;
  }

  //===========================================================================
  // Shutdown MPI
  //===========================================================================

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
