/*
  Generates a Monte Carlo realization of a halo with an embedded
  two-dimensional disk using Jeans' equations.

  Assumptions:

  1) Spherical halo supplied in an input table

  2) Axisymmetric 2-dimensional flat exponential disk

  3) Halo as spherical velocity ellipsoid

  4) Disk as axisymmetric velocity ellipsoid in the plane (that is,
     $\sigma_r = \sigma_\phi$ determined by solving Jeans' equations
     in cylindrical coordinates.

 Updated to use FlatDisk port, Disk2d on 03/23 by MDW

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

#include <config_exp.h>
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
#include <Disk2d.H>
#include <libvars.H>		// Library globals
#include <cxxopts.H>		// Command-line parsing
#include <EXPini.H>		// Ini-style config

#include <norminv.H>

using namespace __EXP__;	// Reference to n-body globals

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
#include <Disk2dHalo.H>
#include <localmpi.H>
#include <Particle.H>

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

  int          NUMR, SCMAP, NUMDF;
  double       RMIN, RCYLMIN, RCYLMAX, SCSPH, RSPHSL, DMFAC, SHFAC, ACYL;
  double       X0, Y0, Z0, U0, V0, W0;
  int          RNUM, PNUM, TNUM, VFLAG, DFLAG;
  bool         expcond, LOGR, CHEBY, SELECT, DUMPCOEF;
  int          CMAPR, CMAPZ, NCHEB, TCHEB, CMTYPE, NDR, NDZ, NHR, NHT, NDP;
  int          NMAXH, NMAXD, NMAXFID, LMAX, MMAX, NUMX, NUMY, NOUT, NODD, DF;
  int          DIVERGE, DIVERGE2, SEED, itmax, nthrds;
  double       DIVERGE_RFAC, DIVERGE_RFAC2;
  double       PPower, R_DF, DR_DF;
  double       Hratio, scale_length, scale_lenfkN;
  double       disk_mass, gas_mass, gscal_length, ToomreQ, Temp, Tmin;
  bool         const_height, images, multi, SVD, DENS, basis, zeropos, zerovel;
  bool         report, ignore, evolved;
  int          nhalo, ndisk, ngparam;
  std::string  hbods, dbods, suffix, centerfile, halofile1, halofile2;
  std::string  cachefile, config, gentype, dtype, dmodel, mtype, ctype;
  
  const std::string mesg("Generates a Monte Carlo realization of a halo with an\n embedded disk using Jeans' equations\n");

  cxxopts::Options options(av[0], mesg);

  options.add_options()
    ("h,help", "Print this help message")
    ("T,template", "Write template options file with current and all default values")
    ("c,config", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("hbods", "The output bodyfile for the halo",
     cxxopts::value<string>(hbods)->default_value("halo.bods"))
    ("dbods", "The output bodyfile for the stellar disc",
     cxxopts::value<string>(dbods)->default_value("disk.bods"))
    ("cachefile", "The cache file for the cylindrical basis",
     cxxopts::value<string>(cachefile)->default_value(".eof_2d_cache"))
    ("ctype", "DiskHalo radial coordinate scaling type (one of: Linear, Log,Rat)",
     cxxopts::value<string>(ctype)->default_value("Log"))
    ("NUMX", "Number of knots in radial dimension of meridional grid",
     cxxopts::value<int>(NUMX)->default_value("256"))
    ("NUMY", "Number of knots in vertical dimension of meridional grid",
     cxxopts::value<int>(NUMY)->default_value("128"))
    ("DIVERGE", "Cusp power-law density extrapolation (0 means off)",
     cxxopts::value<int>(DIVERGE)->default_value("0"))
    ("DIVERGE_RFAC", "Exponent for inner cusp extrapolation",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("1.0"))
    ("DIVERGE2", "Cusp power-law extrapolation for number distribution (0 means off)",
     cxxopts::value<int>(DIVERGE2)->default_value("0"))
    ("DIVERGE_RFAC2", "Exponent for inner cusp extrapolation in number density",
     cxxopts::value<double>(DIVERGE_RFAC2)->default_value("1.0"))
    ("DF", "Use change-over from Jeans to Eddington",
     cxxopts::value<int>(DF)->default_value("1"))
    ("R_DF", "Change-over radius from Jeans to Eddington",
     cxxopts::value<double>(R_DF)->default_value("1.0"))
    ("DR_DF", "Width of change-over from Jeans to Eddington",
     cxxopts::value<double>(DR_DF)->default_value("1.0"))
    ("nhalo", "Number of halo particles",
     cxxopts::value<int>(nhalo)->default_value("100000"))
    ("ndisk", "Number of disk particles",
     cxxopts::value<int>(ndisk)->default_value("100000"))
    ("dtype", "Disk density target type",
     cxxopts::value<string>(dtype)->default_value("exponential"))
    ("LMAX", "Harmonic order for halo expansion",
     cxxopts::value<int>(LMAX)->default_value("18"))
    ("MMAX", "Harmonic order for disk expansion",
     cxxopts::value<int>(MMAX)->default_value("12"))
    ("NOUT", "Number of radial terms in diagnostic basis file for cylinder",
     cxxopts::value<int>(NOUT)->default_value("8"))
    ("NODD", "Number of vertically odd terms in cylindrical expansion",
     cxxopts::value<int>(NODD)->default_value("6"))
    ("NMAXH", "Number of radial terms in spherical expansion",
     cxxopts::value<int>(NMAXH)->default_value("18"))
    ("NMAXD", "Number of radial terms in cylindrical expansion",
     cxxopts::value<int>(NMAXD)->default_value("12"))
    ("NMAXFID", "Number of radial terms in Bessel EOF expansion",
     cxxopts::value<int>(NMAXFID)->default_value("64"))
    ("NUMDF", "Number of knots in Eddington inversion grid",
     cxxopts::value<int>(NUMDF)->default_value("1000"))
    ("VFLAG", "",
     cxxopts::value<int>(VFLAG)->default_value("31"))
    ("DFLAG", "",
     cxxopts::value<int>(DFLAG)->default_value("31"))
    ("expcond", "",
     cxxopts::value<bool>(expcond)->default_value("false"))
    ("report", "",
     cxxopts::value<bool>(report)->default_value("false"))
    ("ignore", "",
     cxxopts::value<bool>(ignore)->default_value("false"))
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
    ("CHEBY", "(boolean) Use Chebyshev smoothing for disc velocities",
     cxxopts::value<bool>(CHEBY)->default_value("false"))
    ("DENS", "(boolean) Output density cylindrical basis tables",
     cxxopts::value<bool>(DENS)->default_value("false"))
    ("zeropos", "(boolean) Zero center of mass",
     cxxopts::value<bool>(zeropos)->default_value("true"))
    ("zerovel", "(boolean) Zero center of velocity",
     cxxopts::value<bool>(zerovel)->default_value("true"))
    ("images", "(boolean) Dump disk basis images",
     cxxopts::value<bool>(images)->default_value("false"))
    ("constheight", "(boolean) Use constant scaleheight for disk",
     cxxopts::value<bool>(const_height)->default_value("true"))
    ("DUMPCOEF", "(boolean) dump coefficients",
     cxxopts::value<bool>(DUMPCOEF)->default_value("false"))
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
    ("disk_mass", "Mass of the stellar disk",
     cxxopts::value<double>(disk_mass)->default_value("0.05"))
    ("scale_length", "Scale length for the realized disk",
     cxxopts::value<double>(scale_length)->default_value("0.01"))
    ("ToomreQ", "Toomre Q parameter for the disk",
     cxxopts::value<double>(ToomreQ)->default_value("1.4"))
    ("RMIN", "Minimum halo radius",
     cxxopts::value<double>(RMIN)->default_value("0.005"))
    ("RCYLMIN", "Minimum disk radius",
     cxxopts::value<double>(RCYLMIN)->default_value("0.001"))
    ("RCYLMAX", "Maximum disk radius, in units of ASCALE",
     cxxopts::value<double>(RCYLMAX)->default_value("20.0"))
    ("SCMAP", "Turn on Spherical SL coordinate mapping (1, 2, 0=off)",
     cxxopts::value<int>(SCMAP)->default_value("1"))
    ("SCSPH", "Scale for Spherical SL coordinate mapping",
     cxxopts::value<double>(SCSPH)->default_value("1.0"))
    ("RSPHSL", "Maximum halo expansion radius",
     cxxopts::value<double>(RSPHSL)->default_value("47.5"))
    ("ACYL", "Radial scale length for disk basis construction",
     cxxopts::value<double>(ACYL)->default_value("0.6"))
    ("DMFAC", "Disk mass scaling factor for spherical deprojection model",
     cxxopts::value<double>(DMFAC)->default_value("1.0"))
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
    ("CMAPZ", "Vertical coordinate mapping type for cylindrical grid  (0=none, 1=rational fct)",
     cxxopts::value<int>(CMAPZ)->default_value("1"))
    ("centerfile", "File containing phase-space center",
     cxxopts::value<std::string>(centerfile))
    ("suffix", "Suffix for output files",
     cxxopts::value<std::string>(suffix)->default_value("diag"))
    ("threads", "Number of threads to run",
     cxxopts::value<int>(nthrds)->default_value("1"))
    ("allow", "Allow multimass algorithm to generature negative masses for testing")
    ("nomono", "Allow non-monotonic mass interpolation")
    ("diskmodel", "Table describing the model for the disk plane")
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
    NOUT = std::min<int>(NOUT, NMAXD);

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
		<< "\t" << "Use parameters read from a YAML config file"  << std::endl
		<< "\t" << av[0] << " --config=gendisk.config"  << std::endl << std::endl
		<< "\t" << "Generate a template YAML config file from current defaults"  << std::endl
		<< "\t" << av[0] << " --template=template.yaml" << std::endl << std::endl;
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
  
  if (vm.count("spline")) {
    SphericalModelTable::linear = 0;
  } else {
    SphericalModelTable::linear = 1;
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
    std::cout << "---- Number of threads=" << numthrd << std::endl
	      << "---- Maximum # threads=" << omp_get_max_threads()
	      << std::endl;
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
  
  int n_particlesH, n_particlesD;
  
  if (suffix.size()>0) {
    hbods = hbods + "." + suffix;
    dbods = dbods + "." + suffix;
  }
  
  // Divvy up the particles by core
  //
  n_particlesH = nhalo/numprocs;
  if (myid==0) n_particlesH = nhalo - n_particlesH*(numprocs-1);
  
  n_particlesD = ndisk/numprocs;
  if (myid==0) n_particlesD = ndisk - n_particlesD*(numprocs-1);
  
#ifdef DEBUG  
  std::cout << "Processor " << myid << ": n_particlesH=" << n_particlesH
	    << std::endl
	    << "Processor " << myid << ": n_particlesD=" << n_particlesD
	    << std::endl
#endif
  
  if (nhalo + ndisk <= 0) {
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
  Disk2dHalo::RDMIN       = std::max<double>(RCYLMIN*scale_length, RMIN);
  Disk2dHalo::RHMIN       = RMIN;
  Disk2dHalo::RHMAX       = RSPHSL;
  Disk2dHalo::RDMAX       = RCYLMAX*scale_length;
  Disk2dHalo::NDR         = NDR;
  Disk2dHalo::NDZ         = NDZ;
  Disk2dHalo::NHR         = NHR;
  Disk2dHalo::NHT         = NHT;
  Disk2dHalo::NDP         = NDP;
  Disk2dHalo::SHFACTOR    = SHFAC;
  Disk2dHalo::COMPRESSION = DMFAC;
  Disk2dHalo::NUMDF       = NUMDF;
  Disk2dHalo::Q           = ToomreQ;
  Disk2dHalo::R_DF        = R_DF;
  Disk2dHalo::DR_DF       = DR_DF;
  Disk2dHalo::SEED        = SEED;
  Disk2dHalo::VFLAG       = static_cast<unsigned int>(DFLAG);
  Disk2dHalo::CHEBY       = CHEBY;
  Disk2dHalo::NCHEB       = NCHEB;

  if (vm.count("itmax"))  Disk2dHalo::ITMAX    = itmax;
  if (vm.count("allow"))  Disk2dHalo::ALLOW    = true;
  if (vm.count("nomono")) Disk2dHalo::use_mono = false;
  if (suffix.size())      Disk2dHalo::RUNTAG   = suffix;
  
  AddDisk::use_mpi      = true;
  AddDisk::Rmin         = RMIN;
  
  //===========================Spherical expansion=============================
  
  // SLGridSph::diverge = DIVERGE;
  // SLGridSph::divergexp = DIVERGE_RFAC;
  
  SphericalSL::RMIN = RMIN;
  SphericalSL::RMAX = RSPHSL;
  SphericalSL::NUMR = NUMR;
  // Create expansion only if needed . . .
  std::shared_ptr<SphericalSL> expandh;
  if (nhalo) {
    expandh = std::make_shared<SphericalSL>(halofile1, nthrds, LMAX, NMAXH, SCMAP, SCSPH);
  }

  //===========================Cylindrical expansion===========================
  //
  // Create YAML db

  YAML::Emitter yml;
  yml << YAML::BeginMap;
  yml << YAML::Key << "acyltbl"   << YAML::Value << ACYL;
  yml << YAML::Key << "rcylmin"   << YAML::Value << RCYLMIN;
  yml << YAML::Key << "rcylmax"   << YAML::Value << RCYLMAX;
  yml << YAML::Key << "scale"     << YAML::Value << scale_length;
  yml << YAML::Key << "numx"      << YAML::Value << NUMX;
  yml << YAML::Key << "numy"      << YAML::Value << NUMY;
  yml << YAML::Key << "numr"      << YAML::Value << NUMR;
  yml << YAML::Key << "cmapR"     << YAML::Value << CMAPR;
  yml << YAML::Key << "cmapZ"     << YAML::Value << CMAPZ;
  yml << YAML::Key << "Mmax"      << YAML::Value << MMAX;
  yml << YAML::Key << "nmax"      << YAML::Value << NMAXD;
  yml << YAML::Key << "nmaxfid"   << YAML::Value << NMAXFID;
  yml << YAML::Key << "logr"      << YAML::Value << LOGR;
  yml << YAML::Key << "cachename" << YAML::Value << cachefile;
  yml << YAML::EndMap;
  
                                // Create expansion only if needed . . .
  std::shared_ptr<Disk2d> expandd;

  if (ndisk)  expandd = std::make_shared<Disk2d>(YAML::Load(yml.c_str()));

  //====================Create the disk & halo model===========================

  std::shared_ptr<Disk2dHalo> diskhalo;

  if (multi) {
    if (myid==0) std::cout << "Initializing for a MULTI-MASS halo . . . " << std::flush;
    diskhalo =
      std::make_shared<Disk2dHalo>
      (expandh, expandd,
       scale_length, disk_mass, ctype,
       halofile1, DIVERGE,  DIVERGE_RFAC,
       halofile2, DIVERGE2, DIVERGE_RFAC2,
       Disk2dHalo::getDiskGenType[gentype]);
    if (myid==0) std::cout << "done" << std::endl;

  } else {

    if (myid==0) std::cout << "Initializing for a SINGLE-MASS halo . . . " << std::flush;
    diskhalo = std::make_shared<Disk2dHalo>
      (expandh, expandd,
       scale_length, disk_mass, ctype, halofile1,
       DF, DIVERGE, DIVERGE_RFAC,
       Disk2dHalo::getDiskGenType[gentype]);
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
  diskhalo->zero_com(zeropos);
  diskhalo->zero_cov(zerovel);
  
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
	if (myid==0) std::cout << "Generating halo phase space . . . " << std::endl;
	diskhalo->set_halo(hparticles, nhalo, n_particlesH);
      } else {
	if (myid==0) std::cout << "Generating halo coordinates . . . " << std::endl;
	diskhalo->set_halo_coordinates(hparticles, nhalo, n_particlesH);
	MPI_Barrier(MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==0) std::cout << "done" << std::endl;
    }
  }

  if (ndisk) {
    if (myid==0) std::cout << "Generating disk coordinates . . . " << std::endl;
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

    // Make dispersion table from particle distribution
    //
    if (evolved) diskhalo->table_halo(hparticles);
  }
  
  if (ndisk) {
    if (myid==0) std::cout << "Beginning disk accumulation . . . " << std::flush;
    expandd->setup_accumulation();
    expandd->accumulate(dparticles);
    MPI_Barrier(MPI_COMM_WORLD);

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

    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      std::cout << "done" << std::endl;
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  //====================Make the phase space velocities========================

  if (!multi and !evolved) {
    if (myid==0) std::cout << "Generating halo velocities . . . " << std::endl;
    diskhalo->set_vel_halo(hparticles);
    if (myid==0) std::cout << "done" << std::endl;
  }
  
  if (myid==0) std::cout << "Generating disk velocities . . . " << std::endl;
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
  diskhalo->profile(outprof, dparticles, 3.0e-3*scale_length, 5.0*scale_length, 100);

  //===========================================================================
  // Shutdown MPI
  //===========================================================================

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
