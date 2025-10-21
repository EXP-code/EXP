/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Use ParticleReader class to compute phase space distribution
 *  difference.  
 *
 *  This version reads phase space from mulitple simulations instances
 *  into a summed paired-difference analysis for improving overall
 *  statistical accuracy.
 *
 *  Orient files are not included (yet) but could easily be added.
 *
 *  Computes surfaces in (E, K) or (I_1, I_2)
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 1/28/22
 *
 ***************************************************************************/

#include <cmath>
#include <cstdlib>

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <Eigen/Eigen>

#include <H5Cpp.h>

#include <ParticleReader.H>
#include <Centering.H>
#include <massmodel.H>
#include <localmpi.H>
#include <EXPini.H>		// Enhanced option parsing
#include <libvars.H>		// EXP library globals
#include <interp.H>
#include "KDE2d.H"		// Kernel density estimation

void p_rec(std::ofstream& out, double E, double K, double V)
{
  out << std::setw(15) << E
      << std::setw(15) << K
      << std::setw(15) << V << std::endl;
}


std::vector<double> p1 {0.0, 0.0, 0.0};	// Phase space #1 center
std::vector<double> p2 {0.0, 0.0, 0.0};	// Phase space #2 center

int
main(int argc, char **argv)
{
  // Kappa offset
  //
  const double ITOL = 1.0e-2;
  const double KTOL = 1.0e-2;

  // MPI initialization
  //
  local_init_mpi(argc, argv);

  // Command-line caching
  //
  std::string cmd_line;
  for (int i=0; i<argc; i++) {
    cmd_line += argv[i];
    cmd_line += " ";
  }

  // Parameter assignment
  //
  double       RMIN;
  double       RMAX;
  double       BMIN;
  double       BMAX;
  double       EMIN;
  double       EMAX;
  double       KMIN;
  double       KMAX;
  double       Sig1, Sig2;
  double       I1min, I1max, I2min, I2max;
  int          NLIMIT;
  int          NSKIP;
  int          NREPORT;
  int          NUMR;
  int          NUM1;
  int          NUM2;
  int          DIVERGE;
  int          nthrds;
  int          WHICHEK;
  int          RNUM;
  int          POSNEG;
  int          TAG;
  int          cstride;
  int          cball=32;
  double       MINBIN;
  bool         SPECIFIC;
  bool         CUMULATE;
  bool         LZDIST;
  bool         LOGMOD;
  bool         meshgrid;
  bool         KDcenter;
  double       DIVERGE_RFAC;
  double       KPOWER;
  std::string  COMP;
  std::string  CURDIR;
  std::string  MODELFILE;
  std::string  OUTFILE;
  std::string  fileType;
  std::string  config;

  std::vector<std::string> INFILE1, INFILE2;

  const char* desc = 
    "=======================================================\n"		\
    "Compute phase-space distributions (DF), phase-space DF \n"         \
    "differences and action changes in phase space for      \n"		\
    "one or more realizations                               \n"		\
    "=======================================================\n"		\
    "   Output file key:\n"						\
    "   ----------------\n"						\
    "   OUTFILE.DN       Counts per bin                 [2d]\n"		\
    "   OUTFILE.DM       Mass per bin                   [2d]\n"		\
    "   OUTFILE.DE       Delta E                        [2d]\n"		\
    "   OUTFILE.DK       Delta L (ang mom mag)          [2d]\n"		\
    "   OUTFILE.DIp      Delta J (phi action)           [2d]\n"         \
    "   OUTFILE.DIr      Delta I (radial action)        [2d]\n"		\
    "   OUTFILE.DJ       Ang mom per bin                [2d]\n"		\
    "   OUTFILE.DKJ      Delta J/J_avg(E)               [2d]\n"		\
    "   OUTFILE.Dm       Mean mass per bin              [2d]\n"		\
    "   OUTFILE.DF       Distribution function          [2d]\n"		\
    "   OUTFILE.ra       Apocentric radius              [2d]\n"         \
    "   OUTFILE.rp       Pericentric radius             [2d]\n"         \
    "   OUTFILE.O1       Radial frequency               [2d]\n"         \
    "   OUTFILE.O2       Azimuthal frequency            [2d]\n"         \
    "   OUTFILE.F1       Two dimensional DF info (1)    [2d]\n"		\
    "   OUTFILE.F2       Two dimensional DF info (2)    [2d]\n"		\
    "   OUTFILE.DR       Run mass, J, Delta J (R)       [1d]\n"		\
    "   OUTFILE.df1      One dimensional DF info (1)    [1d]\n"		\
    "   OUTFILE.df2      One dimensional DF info (2)    [1d]\n"		\
    "   OUTFILE.chk      Orbital element check              \n"		\
    "=======================================================\n"
    " E=Energy, K=J/J_{max}(E)                              \n"
    "=======================================================\n";



  cxxopts::Options options(argv[0], desc);

  options.add_options()
    ("h,help", "This help message")
    ("X,noCommand", "do not save command line")
    ("rmaxf", "Do not extend orbit evaluation beyond the model radius")
    ("relJ", "Compute the relative binned angular momentum")
    ("jaco", "Compute phase-space Jacobian for DF computation")
    ("Emass", "Create energy bins approximately uniform in mass using potential from the mass model")
    ("actions", "Print output in action space rather than E-kappa space.  The default is Energy-Kappa.")
    ("F,filetype", "input file type (one of: PSPout, PSPspl, PSPhdf5, GadgetNative, GadgetHDF5)",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("I1min", "Minimum grid value for E (or I1 for actions)",
     cxxopts::value<double>(I1min))
    ("I1max", "Maximum grid value for E (or I1 for actions)",
     cxxopts::value<double>(I1max))
    ("I2min", "Minimum grid value for K (or I2 for actions)",
     cxxopts::value<double>(I2min))
    ("I2max", "Maximum grid value for K (or I2 for actions)",
     cxxopts::value<double>(I2max))
    ("RMIN", "Minimum model radius",
     cxxopts::value<double>(RMIN)->default_value("0.0"))
    ("RMAX", "Maximum model radius",
     cxxopts::value<double>(RMAX)->default_value("1.0"))
    ("BMIN", "Minimum value in cos(b)",
     cxxopts::value<double>(BMIN)->default_value("-10.0"))
    ("BMAX", "Maximum value in cos(b)",
     cxxopts::value<double>(BMAX)->default_value("10.0"))
    ("EMIN", "Minimum value in phase-space energy",
      cxxopts::value<double>(EMIN)->default_value("-1.0e20"))
    ("EMAX", "Maximum value in phase-space energy",
     cxxopts::value<double>(EMAX)->default_value("1.0e20"))
    ("KMIN", "Minimum value in J/J_max(E)",
     cxxopts::value<double>(KMIN)->default_value("0.0"))
    ("KMAX", "Minimum value in J/J_max(E)",
      cxxopts::value<double>(KMAX)->default_value("1.0"))
    ("NLIMIT", "Number of particles to compare (NLIMIT<0 means all)",
     cxxopts::value<int>(NLIMIT)->default_value("-1"))
    ("NSKIP", "Number of particles to skip before comparing",
     cxxopts::value<int>(NSKIP)->default_value("0"))
    ("CSTRIDE", "Number of particles to stride for center analysis",
     cxxopts::value<int>(cstride)->default_value("10"))
    ("NREPORT", "Interval for reporting processing progress",
     cxxopts::value<int>(NREPORT)->default_value("0"))
    ("NUMR", "Number of radius bins",
     cxxopts::value<int>(NUMR)->default_value("100"))
    ("NUM1", "Number of radial action bins",
     cxxopts::value<int>(NUM1)->default_value("60"))
    ("NUM2", "Number of angular action bins",
     cxxopts::value<int>(NUM2)->default_value("60"))
    ("DIVERGE", "Flag for cusp extrapolation (non-zero means ON)",
     cxxopts::value<int>(DIVERGE)->default_value("0"))
    ("nthrds", "Desired number of OpenMP threads",
     cxxopts::value<int>(nthrds)->default_value("1"))
    ("WHICHEK", "Choose the form of angular momentum indexing phase-space (1=first, 2=second, 3=RMS",
     cxxopts::value<int>(WHICHEK)->default_value("1"))
    ("RNUM", "Size of radial grid for model",
     cxxopts::value<int>(RNUM)->default_value("1000"))
    ("POSNEG", "Use only pos (1) or neg(-1) values for ang. mom. (default=0->all)",
     cxxopts::value<int>(POSNEG)->default_value("0"))
    ("TAG", "Reject particles with integer>0 in attirbute# tag",
     cxxopts::value<int>(TAG)->default_value("-1"))
    ("MINBIN", "Minimum number of particles to bin for output",
     cxxopts::value<double>(MINBIN)->default_value("3.0"))
    ("SPECIFIC", "Normalize by mass per bin",
     cxxopts::value<bool>(SPECIFIC)->default_value("true"))
    ("CUMULATE", "Output cumulative distributions",
     cxxopts::value<bool>(CUMULATE)->default_value("false"))
    ("LZDIST", "Output distribution of L_z in E-K",
     cxxopts::value<bool>(LZDIST)->default_value("false"))
    ("LOGMOD", "Use log(r) scaling for mass model",
     cxxopts::value<bool>(LOGMOD)->default_value("false"))
    ("meshgrid", "Print dimensions for Python meshgrid-like arrays",
     cxxopts::value<bool>(meshgrid)->default_value("true"))
    ("KD", "Use KD to find the density center and use it as the expansion center",
     cxxopts::value<bool>(KDcenter)->default_value("false"))
    ("DIVERGE_RFAC", "Cusp index for extrapolation",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("1.0"))
    ("Kpow", "Create kappa bins with power scaling",
     cxxopts::value<double>(KPOWER)->default_value("1.0"))
    ("Sig1", "Bin smoothing in Dimension 1 (must be > 0)",
     cxxopts::value<double>(Sig1)->default_value("0.0"))
    ("Sig2", "Bin smoothing in Dimension 2 (must be > 0)",
     cxxopts::value<double>(Sig2)->default_value("0.0"))
    ("INFILE1", "Fiducial phase-space file",
     cxxopts::value<std::vector<std::string>>(INFILE1))
    ("INFILE2", "Evolved phase-space file",
     cxxopts::value<std::vector<std::string>>(INFILE2))
    ("CURDIR", "Alternative directory",
     cxxopts::value<string>(CURDIR)->default_value(""))
    ("MODELFILE", "Model file for phase-space orbital values",
     cxxopts::value<string>(MODELFILE)->default_value("SLGridSph.model"))
    ("COMP", "Compute wake for this component name",
     cxxopts::value<std::string>(COMP)->default_value("stars"))
    ("OUTFILE", "Prefix for output files",
     cxxopts::value<std::string>(OUTFILE)->default_value("diffpsP"))
    ("f,input", "Input parameter config file",
     cxxopts::value<std::string>(config))
    ("T,template", "Write template options file with current and all default values")
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // Write YAML template config file and exit
  //
  if (vm.count("template")) {
    // Write template file
    //
    if (myid==0) SaveConfig(vm, options, "template.yaml");
    //
    MPI_Finalize();
    return 0;
  }

  if (vm.count("help")) {
    if (myid==0) std::cout << options.help() << std::endl;
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

  if (myid==0 and vm.count("noCommand")==0) {
    std::string cmdFile = "diffpsp." + OUTFILE + ".cmd_line";
    std::ofstream cmd(cmdFile.c_str());
    if (!cmd) {
      std::cerr << "mssaprof: error opening <" << cmdFile
		<< "> for writing" << std::endl;
    } else {
      cmd << cmd_line << std::endl;
    }
  }

  // Check for that size of file lists match
  //
  if (INFILE1.size() != INFILE2.size()) {
    if (myid==0)
      std::cerr << "Input file vectors must have paired entries"
		<< std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Check for existence of files in INFILE1 and INFILE2 lists
  //
  unsigned bad = 0;
  std::string path1, path2;
  if (fileType == "PSPhdf5") {
    path1 = CURDIR + INFILE1[0];
    std::filesystem::path dir_path = path1;
    if (std::filesystem::is_directory(dir_path)) {
      INFILE1.clear();
      try {
        for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
	  if (std::filesystem::is_regular_file(entry)) {
	    std::string file = entry.path().string();
	    if (H5::H5File::isHdf5(file)) INFILE1.push_back(file);
	  }
        }
      } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "diffpsp error: " << e.what() << std::endl;
      }

      if (INFILE1.size()==0) bad++;
    }
  }
  else {
    for (auto file : INFILE1) {
      if (not std::filesystem::exists(std::filesystem::path(file))) bad++;
    }
  }

  if (bad) {
    if (myid==0)
      std::cerr << "Could not open " << bad
		<< " file(s) in INFILE1 list" << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (fileType == "PSPhdf5") {
    path2 = CURDIR + INFILE2[0];
    std::filesystem::path dir_path = path2;
    if (std::filesystem::is_directory(dir_path)) {
      INFILE2.clear();
      try {
        for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
	  if (std::filesystem::is_regular_file(entry)) {
	    std::string file = entry.path().string();
	    if (H5::H5File::isHdf5(file)) INFILE2.push_back(file);
	  }
        }
      } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "diffpsp error: " << e.what() << std::endl;
      }

      if (INFILE2.size()==0) bad++;
    }
  }
  else {
    for (auto file : INFILE2) {
      if (not std::filesystem::exists(std::filesystem::path(file))) bad++;
    }
  }

  if (bad) {
    if (myid==0)
      std::cerr << "Could not open " << bad
		<< " file(s) in INFILE2 list" << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  bool jaco = false;
  if (vm.count("jaco")) jaco = true;

  bool actions = false;
  if (vm.count("actions")) actions = true;

  // Allocate histograms
  //
  Eigen::MatrixXd histoC, histoM, histoE, histoJ, histoR, histoI, histoT;
  Eigen::MatrixXd histo1, histo2, rapo, rperi, omega1, omega2;

  histoC = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoM = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoE = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoJ = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoI = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoR = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoT = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histo1 = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histo2 = Eigen::MatrixXd::Zero(NUM1, NUM2);
  rapo   = Eigen::MatrixXd::Zero(NUM1, NUM2);
  rperi  = Eigen::MatrixXd::Zero(NUM1, NUM2);
  omega1 = Eigen::MatrixXd::Zero(NUM1, NUM2);
  omega2 = Eigen::MatrixXd::Zero(NUM1, NUM2);

  Eigen::VectorXd histoP, histoL, histPr, histLr, histoS, histoN;

  histoP = Eigen::VectorXd::Zero(NUMR);
  histoL = Eigen::VectorXd::Zero(NUMR);
  histPr = Eigen::VectorXd::Zero(NUMR);
  histLr = Eigen::VectorXd::Zero(NUMR);
  histoS = Eigen::VectorXd::Zero(NUMR);
  histoN = Eigen::VectorXd::Zero(NUMR);

  // One-d histograms in energy or J depending on coordinates
  //
  std::vector<Eigen::VectorXd> histo1_1d(2), histo2_1d(2);

  histo1_1d[0] = Eigen::VectorXd::Zero(NUM1);
  histo2_1d[0] = Eigen::VectorXd::Zero(NUM1);
  histo1_1d[1] = Eigen::VectorXd::Zero(NUM2);
  histo2_1d[1] = Eigen::VectorXd::Zero(NUM2);

  double dR, rhmin, rhmax;
  bool LOGR;
  if (RMIN>1.0e-08) {
    LOGR = true;
    rhmin = log(RMIN);
    rhmax = log(RMAX);
  } else {
    LOGR = false;
    rhmin = RMIN;
    rhmax = RMAX;
  }
  dR = (rhmax - rhmin)/(NUMR-1);

  //
  // Open output file
  //
    
  const int nfiles = 21;
  const char *suffix[] = {
    ".DM", 			// Mass per bin (E, K)		#0
    ".DE", 			// Delta E(E, K)		#1
    ".DK", 			// Delta J(E, K)		#2
    ".DJ", 			// Ang mom per bin (E, K)	#3
    ".DIr", 			// Rad. act. per bin (E, K)	#4
    ".DIp", 			// Phi. act. per bin (E, K)	#5
    ".DKJ", 			// Delta J(E, K)/J_avg(E)	#6
    ".Dm", 			// Mean mass per bin (E, K)	#7
    ".DF", 			// Delta DF (E, K)              #8
    ".ra", 			// Apocentric radius (E, K)     #9
    ".rp", 			// Pericentric radius (E, K)    #10
    ".O1", 			// Radial frequency (E, K)      #11
    ".O2", 			// Azimuthal frequency (E, K)   #12
    ".F1", 			// DF for PS 1 (E, K)           #13
    ".F2", 			// DF for PS 2 (E, K)           #14
    ".Df", 			// Delta DF/F (E, K)            #15
    ".DN", 			// Counts per (E, K) bin        #16
    ".DR", 			// Run mass, J, Delta J (R)	#17
    ".chk",			// Orbital element check	#18
    ".df1",			// One dimensional DF info	#19
    ".df2"			// One dimensional DF info	#20
  };
  std::vector<string> filename(nfiles);
  for (int i=0; i<nfiles; i++) filename[i] = OUTFILE + suffix[i];

  std::vector<std::ofstream> out(nfiles);
  
  bool good = true;
  if (myid==0) {
    for (int i=0; i<nfiles; i++) {
    
      out[i].open(filename[i].c_str());

      if (!out[i]) {
	std::cerr << "Error opening <" << filename[i] << ">" << std::endl;
	good = false;
      }
    }
  }
  
  MPI_Bcast(&good, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

  if (not good) {
    MPI_Finalize();
    exit(-1);
  }

  //============================================================
  // Build model
  //============================================================

  PR::PRptr psp;

  auto hmodel = std::make_shared<SphericalModelTable>(MODELFILE, DIVERGE, DIVERGE_RFAC);

  if (vm.count("rmaxf")) {
    SphericalOrbit::RMAXF = 1.0;
    std::cout << "RMAXF=" << SphericalOrbit::RMAXF << std::endl;
  }

  SphericalOrbit orb(hmodel);

  double Emin = std::max<double>(hmodel->get_pot(hmodel->get_min_radius()), EMIN);
  double Emax = std::min<double>(hmodel->get_pot(hmodel->get_max_radius()), EMAX);

  std::vector<double> Ebins;
  if (vm.count("Emass")) {
    // Make bins
    const int Ntmp = 1000;
    std::vector<double> mvec(Ntmp), pvec(Ntmp);
    double rmin = hmodel->get_min_radius();
    double rmax = hmodel->get_max_radius();
    double dr   = (rmax - rmin)/(Ntmp-1);
    for (int i=0; i<Ntmp; i++) {
      double r = rmin + dr*i;
      mvec[i] = hmodel->get_mass(r);
      pvec[i] = hmodel->get_pot(r);
    }
    // Get mass range for energy range
    Linear1d ENGY(pvec, mvec);
    Linear1d MASS(mvec, pvec);

    Ebins.resize(NUM1+1);
    double Mmin = ENGY.eval(Emin), Mmax = ENGY.eval(Emax);
    double dM = (Mmax - Mmin)/NUM1;
    for (int i=0; i<=NUM1; i++) Ebins[i] = MASS.eval(Mmin+dM*i);
  }

  if (actions) {

    Emin *= 1.0 - KTOL;
    Emax *= 1.0 + KTOL;

    orb.new_orbit(Emin, 1.0 - KTOL);
    if (vm.count("I1min")==0) I1min = orb.get_action(0);

    orb.new_orbit(Emin, KTOL);
    if (vm.count("I2min")==0) I2min = orb.get_action(1);

    orb.new_orbit(Emax, KTOL);
    if (vm.count("I1max")==0) I1max = orb.get_action(0);

    orb.new_orbit(Emax, 1.0 - KTOL);
    if (vm.count("I2max")==0) I2max = orb.get_action(1);

    if (myid==0)
      std::cout << std::endl
		<< std::string(40, '-')  << std::endl
		<< "---- Range limits"  << std::endl
		<< std::string(40, '-')  << std::endl
		<< "-- I1min: " << I1min << std::endl
		<< "-- I1max: " << I1max << std::endl
		<< "-- I2min: " << I2min << std::endl
		<< "-- I2max: " << I2max << std::endl
		<< std::string(40, '-')  << std::endl;
  } else {
    if (vm.count("I1min")==0) I1min = Emin;
    if (vm.count("I1max")==0) I1max = Emax;
    if (vm.count("I2min")==0) I2min = std::max<double>(KTOL, KMIN);
    if (vm.count("I2max")==0) I2max = std::min<double>(1.0 - KTOL, KMAX);

    if (myid==0)
      std::cout << std::endl
		<< std::string(40, '-')  << std::endl
		<< "---- Range limits"  << std::endl
		<< std::string(40, '-')  << std::endl
		<< "-- Emin: " << I1min << std::endl
		<< "-- Emax: " << I1max << std::endl
		<< "-- Kmin: " << I2min << std::endl
		<< "-- Kmax: " << I2max << std::endl
		<< std::string(40, '-')  << std::endl;
  }

  double d1 = (I1max - I1min) / NUM1;
  double d2 = (I2max - I2min) / NUM2;
  double Ir1, Ip1, Ir2, Ip2, I1, I2;

  double d3 = d2;
  bool Kpow = false;
  if (not actions and vm.count("Kpow")) {
    Kpow = true;
    d3 = (pow(I2max, KPOWER) - pow(I2min, KPOWER)) / NUM2;
  }


  //=======================
  // Open distribution file
  //=======================

  std::ofstream dout;
  if (LZDIST and myid==0) {
    string DSTFILE = OUTFILE + ".dist";
    dout.open(DSTFILE.c_str());
    if (!dout) {
      cerr << "Couldn't open <" << DSTFILE << ">" << std::endl;
      LZDIST = false;
    }
  }


  // Diagnostic values
  //
  int reject=0, N=0, total=0, rover=0, emiss=0, pmiss=0, Ntot=0;

  // Times for the PSP snaps
  //
  double initl_time, final_time;

  // Number of paths
  //
  int npath1 = 1;
  if (fileType != "PSPhdf5") {
    npath1 = INFILE1.size();
  }

  // Iterate through file list
  //
  for (size_t n=0; n<npath1; n++) {

    PR::PRptr psp1, psp2;

    try {

      if (fileType == "PSPhdf5") {
	psp1 = PR::ParticleReader::createReader(fileType, INFILE1, myid, true);
      } else {
	psp1 = PR::ParticleReader::createReader(fileType, {INFILE1[n]}, myid, true);
      }
  
      initl_time = psp1->CurrentTime();

      psp1->SelectType(COMP);

      if (myid==0) {
	std::cout << std::endl << std::string(40, '-') << std::endl;
	if (fileType == "PSPhdf5")
	  std::cout << "Path 1: " << path1 << std::endl;
	else
	  std::cout << "File 1: " << CURDIR + INFILE1[n] << std::endl;
	std::cout << "Found dump at time: " << initl_time << std::endl;
      }
    }
    catch (const std::runtime_error& error) {
      if (myid==0)
	std::cerr << "diffpsp: error opening snapshot in file <"
		  << INFILE1[n] << ">" << std::endl
		<< "diffpsp: " << error.what()
		<< std::endl;
      MPI_Finalize();
      exit(-1);
    }
    
    try {
      if (fileType == "PSPhdf5") {
	std::string path = CURDIR + INFILE2[0];
	psp2 = PR::ParticleReader::createReader(fileType, INFILE2, myid, true);
      } else {
	psp2 = PR::ParticleReader::createReader(fileType, {INFILE2[n]}, myid, true);
      }
  
      final_time = psp2->CurrentTime();

      psp2->SelectType(COMP);

      if (myid==0) {
	if (fileType == "PSPhdf5")
	  std::cout << "Path 2: " << path2 << endl;
	else
	  std::cout << "File 2: " << INFILE2[n] << endl;
	std::cout << "Found dump at time: " << final_time << std::endl;
	std::cout << std::string(40, '-') << std::endl;
      }
    }
    catch (const std::runtime_error& error) {
      if (myid==0)
      std::cerr << "diffpsp: error opening snapshot in file <"
		<< INFILE2[n] << ">"<< std::endl
		<< "diffpsp: " << error.what()
		<< std::endl;
      MPI_Finalize();
      exit(-1);
    }
    

    //============================================================
    // Default centers using KD?
    //============================================================
    
    if (KDcenter) {

      if (myid==0) {		// stdout info
	std::cout << std::endl
		  << std::string(52, '-')     << std::endl
		  << "---- Center estimation" << std::endl
		  << std::string(52, '-')     << std::endl;
      }

      p1 = Utility::getDensityCenter(psp1, cstride, 0, cball);
      p2 = Utility::getDensityCenter(psp2, cstride, 0, cball);

      if (myid==0) {		// stdout info
	std::cout << "#:  "
		  << std::setw(16) << std::right << "X    "
		  << std::setw(16) << std::right << "Y    "
		  << std::setw(16) << std::right << "Z    "
		  << std::endl
		  << "--  " << std::setfill('-')
		  << std::setw(16) << std::left << "  "
		  << std::setw(16) << std::left << "  "
		  << std::setw(16) << std::left << "  "
		  << std::endl << std::setfill(' ');
	std::cout << "0:  ";
	for (auto v : p1) std::cout << std::setw(16) << std::right << v;
	std::cout << std::endl << "1:  ";
	for (auto v : p2) std::cout << std::setw(16) << std::right << v;
	std::cout << std::endl << std::string(52, '-') << std::endl
		  << std::endl;
      }
    }

    //============================================================
    // Do difference
    //============================================================
    
    struct SPS
    {
      double pos[3], vel[3];
    } tps;
    
    std::map<int, SPS> ph;

    int N;

    {
      std::vector<double> pos0, vel0, pos1, vel1;
      std::vector<int> ind0, ind1;

      if (myid==0 and NREPORT) {
	std::cout << "Initial map creation" << std::endl;
      }

      N = 0;
      for (auto pp=psp1->firstParticle(); pp!=0; pp=psp1->nextParticle()) {
	for (int k=0; k<3; k++) {
	  // Add particles to map
	  tps.pos[k] = pp->pos[k];
	  tps.vel[k] = pp->vel[k];
	  // Pack particles for MPI
	  if (numprocs>1) {
	    pos0.push_back(pp->pos[k]);
	    vel0.push_back(pp->vel[k]);
	  }
	}
	ph[pp->indx] = tps;
	if (numprocs>1) ind0.push_back(pp->indx);

	if (myid==0 and NREPORT) {
	  if (!((N+1)%NREPORT))
	    std::cout << "\rProcessed: " 
		      << std::setw(10) << (N+1)*numprocs << std::flush;
	  N++;
	}
      }

      if (myid==0 and NREPORT) {
	std::cout << std::endl << std::endl
		  << "Particle map exchange..." << std::endl;
      }

      // Exchange particles with other processes
      //
      if (numprocs>1) {
	for (int n=0; n<numprocs; n++) {
	  if (n==myid) {
	    N = ind0.size();
	    MPI_Bcast(&N, 1, MPI_INT, n, MPI_COMM_WORLD);
	    MPI_Bcast(ind0.data(),   N, MPI_INT,    n, MPI_COMM_WORLD);
	    MPI_Bcast(pos0.data(), 3*N, MPI_DOUBLE, n, MPI_COMM_WORLD);
	    MPI_Bcast(vel0.data(), 3*N, MPI_DOUBLE, n, MPI_COMM_WORLD);
	  } else {
	    MPI_Bcast(&N, 1, MPI_INT, n, MPI_COMM_WORLD);
	    ind1.resize(N);
	    pos1.resize(3*N);
	    vel1.resize(3*N);
	    MPI_Bcast(ind1.data(),   N, MPI_INT,    n, MPI_COMM_WORLD);
	    MPI_Bcast(pos1.data(), 3*N, MPI_DOUBLE, n, MPI_COMM_WORLD);
	    MPI_Bcast(vel1.data(), 3*N, MPI_DOUBLE, n, MPI_COMM_WORLD);

	    // Load the map
	    for (int i=0; i<N; i++) {
	      for (int k=0; k<3; k++) {
		tps.pos[k] = pos1[i*3+k];
		tps.vel[k] = vel1[i*3+k];
	      }
	      ph[ind1[i]] = tps;
	    }
	  }

	  if (myid==0 and NREPORT) {
	    std::cout << "\rNode " << std::setw(4) << n
		      << ": " << N << std::flush;
	  }

	}
	// END: process loop

	if (myid==0 and NREPORT and numprocs>1) std::cout << std::endl;
      }
      // END: MPI exchange
    }

    // Begin particle difference loop
    //
    if (myid==0 and NREPORT) {
      std::cout << std::endl << "Particle differencing..." << std::endl;
    }

    N = 0;
    for (auto pp=psp2->firstParticle(); pp!=0; pp=psp2->nextParticle()) {
    
      auto ip = ph.find(pp->indx);
    
      if (ip != ph.end()) {
      
	double angmom1[3], angmom2[3];
	double p10[3], p20[3], v10[3], v20[3];
      
	if (TAG>=0 && pp->iattrib.size()>TAG) {
	  if (pp->iattrib[TAG]>0) continue;
	}
	
	for (int k=0; k<3; k++) p10[k] = ip->second.pos[k] - p1[k];
	for (int k=0; k<3; k++) p20[k] = pp->pos[k] - p2[k];
	for (int k=0; k<3; k++) v10[k] = ip->second.vel[k];
	for (int k=0; k<3; k++) v20[k] = pp->vel[k];
      
	double rr1=0, vv1=0, rr2=0, vv2=0;
	for (int k=0; k<3; k++) {
	  rr1 += p10[k]*p10[k];
	  vv1 += v10[k]*v10[k];
	  rr2 += p20[k]*p20[k];
	  vv2 += v20[k]*v20[k];
	}
      
	if (rr1 <= RMAX*RMAX && rr1 >= RMIN*RMIN) {
	
	  double E1 = 0.5*vv1 + hmodel->get_pot(sqrt(rr1));

	  if (E1 < Emin or E1 > Emax) {emiss++; continue;}

	  double E2 = 0.5*vv2 + hmodel->get_pot(sqrt(rr2));
	  
	  if (E2 < Emin or E2 > Emax) {emiss++; continue;}
	  
	  angmom1[0] = p10[1]*v10[2] - p10[2]*v10[1];
	  angmom1[1] = p10[2]*v10[0] - p10[0]*v10[2];
	  angmom1[2] = p10[0]*v10[1] - p10[1]*v10[0];
	
	  double cosb = angmom1[2]/sqrt(angmom1[0]*angmom1[0] +
					angmom1[1]*angmom1[1] +
					angmom1[2]*angmom1[2] );
	
	  if (POSNEG>0 && angmom1[2]<0.0 ||
	      POSNEG<0 && angmom1[2]>0.0) continue;

	  if (cosb<BMIN || cosb>BMAX) continue;
	
	  angmom2[0] = p20[1]*v20[2] - p20[2]*v20[1];
	  angmom2[1] = p20[2]*v20[0] - p20[0]*v20[2];
	  angmom2[2] = p20[0]*v20[1] - p20[1]*v20[0];
	  
	  double jj = 0.0, dj = 0.0, j1 = 0.0, j2 = 0.0;
	  for (int k=0; k<3; k++) {
	    if (WHICHEK==1)
	      jj += angmom1[k]*angmom1[k];
	    else if (WHICHEK==2)
	      jj += angmom2[k]*angmom2[k];
	    else
	      jj += 0.5*(angmom1[k]*angmom1[k] + angmom2[k]*angmom2[k]);
	  
	    j1 += angmom1[k]*angmom1[k];
	    j2 += angmom2[k]*angmom2[k];
	    dj += (angmom1[k] - angmom2[k])*(angmom1[k] - angmom2[k]);
	  }
	
	  try {
	    orb.new_orbit(E1, 0.5);
	    double K1 = sqrt(j1)/orb.Jmax();

	    if (K1>1.0-KTOL) throw std::runtime_error("K1 > 1-KTOL");
	    if (K1<KTOL)     throw std::runtime_error("K1 < KTOL");
	    
	    orb.new_orbit(E1, K1);
	  
	    Ir1 = orb.get_action(0);
	    Ip1 = orb.get_action(1);
	  
	    orb.new_orbit(E2, 0.5);
	    double K2 = sqrt(j2)/orb.Jmax();

	    if (K2>1.0-KTOL) throw std::runtime_error("K2 > 1-KTOL");
	    if (K2<KTOL)     throw std::runtime_error("K2 < KTOL");
	  
	    orb.new_orbit(E2, K2);

	    double Ir2 = orb.get_action(0);
	    double Ip2 = orb.get_action(1);
	  
	    if (myid==0) {
	      if (WHICHEK & 1) {
		if (K1>1.0 || K1<0.0)
		  out[12] << setw(15) << E1 << setw(15) << K1
			  << setw(15) << E2 << setw(15) << sqrt(j1)
			  << setw(15) << sqrt(rr1) << setw(15) << sqrt(rr2)
			  << setw(5) << 1 << endl;
	      }
	    
	      if (WHICHEK & 2) {
		if (K2>1.0 || K2<0.0)
		  out[12] << setw(15) << E2 << setw(15) << K2
			  << setw(15) << E1 << setw(15) << sqrt(j2)
			  << setw(15) << sqrt(rr1) << setw(15) << sqrt(rr2)
			  << setw(5) << 2 << endl;
	      }
	    }
	  
	    double EE, KK;
	    if (WHICHEK == 1) {
	      EE = E1;
	      KK = K1;
	    }
	    else if (WHICHEK == 2) {
	      EE = E2;
	      KK = K2;
	    }
	    else {
	      EE = 0.5*(E1 + E2);
	      KK = 0.5*(K1 + K2);
	    }
	    
	    if (EE > Emax or EE < Emin) continue;
	    if (KK > 1.0 - KTOL or KK < KTOL) continue;

	    double I1_1, I2_1, I1_2, I2_2, ra, rp, O1, O2;
	    int i1, i2, i11, i12, i21, i22;

	    orb.new_orbit(EE, KK);
	    ra = orb.apo();
	    rp = orb.peri();
	    O1 = orb.get_freq(0);
	    O2 = orb.get_freq(1);

	    if (actions) {

	      // Get the actions
	      //
	      I1 = orb.get_action(0);
	      I2 = orb.get_action(1);

	      if (WHICHEK == 1) {
		I1_1 = I1;
		I2_1 = I2;
		orb.new_orbit(E2, K2);
		I1_2 = orb.get_action(0);
		I2_2 = orb.get_action(1);
	      } else if (WHICHEK == 2) {
		I1_2 = I1;
		I2_2 = I2;
		orb.new_orbit(E1, K1);
		I1_1 = orb.get_action(0);
		I2_1 = orb.get_action(1);
	      } else {
		orb.new_orbit(E1, K1);
		I1_1 = orb.get_action(0);
		I2_1 = orb.get_action(1);
		orb.new_orbit(E2, K2);
		I1_2 = orb.get_action(0);
		I2_2 = orb.get_action(1);
	      }

	      if (I1_1 < I1min) throw std::runtime_error("I1 < I1min [1]");
	      if (I1_1 > I1max) throw std::runtime_error("I1 > I1max [1]");
	      if (I2_1 < I2min) throw std::runtime_error("I2 < I2min [1]");
	      if (I2_1 > I2max) throw std::runtime_error("I2 > I2max [1]");

	      if (I1_2 < I1min) throw std::runtime_error("I1 < I1min [2]");
	      if (I1_2 > I1max) throw std::runtime_error("I1 > I1max [2]");
	      if (I2_2 < I2min) throw std::runtime_error("I2 < I2min [2]");
	      if (I2_2 > I2max) throw std::runtime_error("I2 > I2max [2]");

	      i1 = (int)floor( (I1 - I1min) / d1 );
	      i1 = std::max<int>(i1, 0);
	      i1 = std::min<int>(i1, NUM1-1);
	    
	      i2 = (int)floor( (I2 - I2min) / d2 );
	      i2 = std::max<int>(i2, 0);
	      i2 = std::min<int>(i2, NUM2-1);
	    
	      if (WHICHEK == 1) {
		i11 = i1;
		i21 = i2;

		i12 = (int)floor( (I1_2 - I1min) / d1 );
		i12 = std::max<int>(i12, 0);
		i12 = std::min<int>(i12, NUM1-1);
		
		i22 = (int)floor( (I2_2 - I2min) / d2 );
		i22 = std::max<int>(i22, 0);
		i22 = std::min<int>(i22, NUM2-1);
	      } else if (WHICHEK == 2) {
		i12 = i1;
		i22 = i2;

		i11 = (int)floor( (I1_1 - I1min) / d1 );
		i11 = std::max<int>(i11, 0);
		i11 = std::min<int>(i11, NUM1-1);
		
		i21 = (int)floor( (I2_1 - I2min) / d2 );
		i21 = std::max<int>(i21, 0);
		i21 = std::min<int>(i21, NUM2-1);
	      } else {
		i11 = (int)floor( (I1_1 - I1min) / d1 );
		i11 = std::max<int>(i11, 0);
		i11 = std::min<int>(i11, NUM1-1);
	    
		i21 = (int)floor( (I2_1 - I2min) / d2 );
		i21 = std::max<int>(i21, 0);
		i21 = std::min<int>(i21, NUM2-1);
		
		i12 = (int)floor( (I1_2 - I1min) / d1 );
		i12 = std::max<int>(i12, 0);
		i12 = std::min<int>(i12, NUM1-1);
	    
		i22 = (int)floor( (I2_2 - I2min) / d2 );
		i22 = std::max<int>(i22, 0);
		i22 = std::min<int>(i22, NUM2-1);
	      }
	      
	    } else {

	      if (E1 < Emin) throw std::runtime_error("E1 < Emin");
	      if (E1 > Emax) throw std::runtime_error("E1 > Emax");
	      if (E2 < Emin) throw std::runtime_error("E2 < Emin");
	      if (E2 > Emax) throw std::runtime_error("E2 > Emax");

	      if (K1 < KTOL) throw std::runtime_error("K1 < KTOL");
	      if (K2 < KTOL) throw std::runtime_error("K2 < KTOL");

	      if (K1 > 1.0 - KTOL) throw std::runtime_error("K1 > 1-KTOL");
	      if (K2 > 1.0 - KTOL) throw std::runtime_error("K2 > 1-KTOL");

	      i11 = (int)floor( (E1 - Emin) / d1 );
	      if (i11 < 0)     throw std::runtime_error("i11 < 0");
	      if (i11 >= NUM1) throw std::runtime_error("i11 >= NUM1");
	    
	      i21 = (int)floor( K1 / d2 );
	      if (i21 < 0)     throw std::runtime_error("i21 < 0");
	      if (i21 >= NUM2) throw std::runtime_error("i21 >= NUM2");

	      i12 = (int)floor( (E2 - Emin) / d1 );
	      if (i12 < 0)     throw std::runtime_error("i12 < 0");
	      if (i12 >= NUM1) throw std::runtime_error("i12 >= NUM1");
	  
	      i22 = (int)floor( K2 / d2 );
	      if (i22 < 0)     throw std::runtime_error("i22 < 0");
	      if (i22 >= NUM2) throw std::runtime_error("i22 >= NUM2");

	      if (Ebins.size()) {
		auto it = std::lower_bound(Ebins.begin(), Ebins.end(), EE);
		if (it != Ebins.begin()) it--;
		i1 = it - Ebins.begin();
		i1 = max<int>(i1, 0);
		i1 = min<int>(i1, NUM1-1);
	      } else {
		i1 = (int)floor( (EE - Emin) / d1 );
		i1 = max<int>(i1, 0);
		i1 = min<int>(i1, NUM1-1);
	      }
	    
	      if (Kpow) {
		i2 = (int)floor( pow(KK, KPOWER) / d3 );
		i2 = max<int>(i2, 0);
		i2 = min<int>(i2, NUM2-1);
	      } else {
		i2 = (int)floor( KK / d2 );
		i2 = max<int>(i2, 0);
		i2 = min<int>(i2, NUM2-1);
	      }
	    }
	    
	    double L1 = 0.0, L2 = 0.0;
	    for (int k=0; k<3; k++) {
	      L1 += angmom1[k]*angmom1[k];
	      L2 += angmom2[k]*angmom2[k];
	    }
	    L1 = sqrt(L1);
	    L2 = sqrt(L2);
	  
	    histoC(i1, i2) += 1;
	    histoM(i1, i2) += pp->mass;
	    histoE(i1, i2) += pp->mass*(E1 - E2);
	    histoJ(i1, i2) += pp->mass*(L1 - L2);
	    histoR(i1, i2) += pp->mass*(Ir2 - Ir1);
	    histoI(i1, i2) += pp->mass*(Ip2 - Ip1);
	    histoT(i1, i2) += pp->mass*L1;
	    
	    rapo  (i1, i2) += pp->mass*ra;
	    rperi (i1, i2) += pp->mass*rp;

	    omega1(i1, i2) += pp->mass*O1;
	    omega2(i1, i2) += pp->mass*O2;

	    histo1(i11, i21) += pp->mass;
	    histo2(i12, i22) += pp->mass;
	    
	    histo1_1d[0](i11) += pp->mass;
	    histo2_1d[0](i12) += pp->mass;
	    histo1_1d[1](i21) += pp->mass;
	    histo2_1d[1](i22) += pp->mass;
	      
	    if (LZDIST and myid==0) {
	      if (KK>KMIN && KK<KMAX && EE>EMIN && EE<EMAX)
		dout << setw(15) << EE
		     << setw(15) << KK
		     << setw(15) << angmom1[2]
		     << setw(15) << angmom2[2]
		     << setw(15) << angmom2[2] - angmom1[2]
		     << endl;
	    }
	  
	    int ir;
	    if (LOGR) 
	      ir = (int)floor( (log(sqrt(rr1)) - rhmin) / dR );
	    else
	      ir = (int)floor( (sqrt(rr1) - rhmin) / dR );
	    ir = max<int>(ir, 0);
	    ir = min<int>(ir, NUMR-1);
	  
	    histoP[ir] += pp->mass*(E1 - E2);
	    histoL[ir] += pp->mass*(angmom2[2] - angmom1[2]);
	    histPr[ir] += pp->mass*(E1 - E2)*2.0/sqrt(E1*E1 + E2*E2);
	    histLr[ir] += pp->mass*(angmom2[2] - angmom1[2])*2.0/
	      sqrt(angmom1[2]*angmom1[2] + angmom2[2]*angmom2[2]);
	    histoN[ir] += pp->mass;
	    histoS[ir] += pp->mass*L1;
	    
	    total++;
	  }
	  catch (const std::runtime_error& error) {
	    if (false)		// Verbose output
	       std::cout << "error [" << reject << "]: "
			 << error.what() << std::endl;
	    reject++;		// Tally grid rejections
	  }
	  Ntot++;
	} else rover++;

      } else pmiss++;
      
      if (myid==0 and NREPORT) {
	if (!((N+1)%NREPORT))
	  std::cout << "\rProcessed: " 
		    << std::setw(10) << (N+1)*numprocs << std::flush;
	N++;
      }

    }
    // END: particle loop

    if (myid==0 and NREPORT)
      std::cout << std::endl << std::string(40, '-') << std::endl;
  }
  // END: file loop
  
  
  // Send to root
  //
  if (myid) {
    MPI_Reduce(&total,        0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rover,        0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&pmiss,        0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&emiss,        0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Ntot,         0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&reject,       0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoC.data(), 0, histoC.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoM.data(), 0, histoM.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rapo.data(),   0, rapo.size(),   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rperi.data(),  0, rperi.size(),  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(omega1.data(), 0, omega1.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(omega2.data(), 0, omega2.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo1.data(), 0, histo1.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo2.data(), 0, histo2.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoE.data(), 0, histoE.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoJ.data(), 0, histoJ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoR.data(), 0, histoR.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoI.data(), 0, histoI.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoT.data(), 0, histoT.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoP.data(), 0, histoP.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoL.data(), 0, histoL.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histPr.data(), 0, histPr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histLr.data(), 0, histLr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoS.data(), 0, histoS.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoN.data(), 0, histoN.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo1_1d[0].data(), 0, histo1_1d[0].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo2_1d[0].data(), 0, histo2_1d[0].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo1_1d[1].data(), 0, histo1_1d[1].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo2_1d[1].data(), 0, histo2_1d[1].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  // Receive at root
  //
  else {
    MPI_Reduce(MPI_IN_PLACE, &total,  1,                   MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &rover,  1,                   MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &pmiss,  1,                   MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &emiss,  1,                   MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &Ntot,   1,                   MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &reject, 1,                   MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoC.data(), histoC.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoM.data(), histoM.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, rapo.data(),   rapo.size(),   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, rperi.data(),  rperi.size(),  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, omega1.data(), omega1.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, omega2.data(), omega2.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo1.data(), histo1.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo2.data(), histo2.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoE.data(), histoE.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoJ.data(), histoJ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoR.data(), histoR.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoI.data(), histoI.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoT.data(), histoT.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoP.data(), histoP.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoL.data(), histoL.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histPr.data(), histPr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histLr.data(), histLr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoS.data(), histoS.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoN.data(), histoN.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo1_1d[0].data(), histo1_1d[0].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo2_1d[0].data(), histo2_1d[0].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo1_1d[1].data(), histo1_1d[1].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo2_1d[1].data(), histo2_1d[1].size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    std::cout << std::endl
	      << std::setw(40) << std::left
	      << "Total number of particles processed: " << total << std::endl
	      << std::setw(40) << std::left
	      << "Out of radial bounds: "                << rover << std::endl
	      << std::setw(40) << std::left
	      << "Out of energy bounds: "                << emiss << std::endl
	      << std::setw(40) << std::left
	      << "Index cache miss: "                    << pmiss << std::endl;

    if (reject) std::cout << std::setw(40) << std::left
			  << "Orbit off-grid failures: " << reject
			  << "/" << Ntot << " states" << std::endl << std::endl;
    std::cout << std::endl;

    
    Eigen::VectorXd I1avg  = Eigen::VectorXd::Zero(NUM1);
    Eigen::VectorXd I2avg  = Eigen::VectorXd::Zero(NUM2);
    Eigen::VectorXd I1mass = Eigen::VectorXd::Zero(NUM1);
    Eigen::VectorXd I2mass = Eigen::VectorXd::Zero(NUM2);

    // Delta DF
    Eigen::MatrixXd histoF  = histo2 - histo1;
    Eigen::MatrixXd histoDF = histoF;

    std::vector<Eigen::VectorXd> histoF_1d =
      {histo2_1d[0] - histo1_1d[0], histo2_1d[1] - histo1_1d[1]};
    std::vector<Eigen::VectorXd> histoDF_1d = histoF_1d;
    
    // Relative Delta DF
    for (size_t j=0; j<histoF.cols(); j++) { // Footron order...
      for (size_t i=0; i<histoF.rows(); i++) {
	if (histo1(i, j)>0.0) histoDF(i, j) = histoF(i, j)/histo1(i, j);
	else                  histoDF(i, j) = 0.0;
      }
    }
	
    for (int i=0; i<2; i++) {
      for (size_t j=0; j<histoF_1d[i].size(); j++) { // 1d array
	if (histo1_1d[i](j)>0.0)
	  histoDF_1d[i](j) = histoF_1d[i](j)/histo1_1d[i](j);
	else
	  histoDF_1d[i](j) = 0.0;
      }
    }

    double totMass = 0.0;
    
    for (int i=0; i<NUM1; i++) {
      for (int j=0; j<NUM2; j++) {
	I1avg[i]  += histoT(i, j);
	I1mass[i] += histoM(i, j);
	totMass   += histoM(i, j);
      }
      if (I1mass[i] > 0.0) I1avg[i] /= I1mass[i];
      else                 I1avg[i]  = 0.0;
    }
    
    for (int j=0; j<NUM2; j++) {
      for (int i=0; i<NUM1; i++) {
	I2avg[j]  += histoT(i, j);
	I2mass[j] += histoM(i, j);
      }
      if (I2mass[j] > 0.0) I2avg[j] /= I2mass[j];
      else                 I2avg[j]  = 0.0;
    }
    
    double mfac = d1*d2;
    
    if (meshgrid)
      for (int k=0; k<17; k++) out[k] << std::setw(8) << NUM1
				      << std::setw(8) << NUM2
				      << std::endl;
    bool relJ = false;
    if (vm.count("relJ")) relJ = true;

    // Smooth arrays
    //
    if (Sig1 > 0.0 and Sig2 > 0.0) {

      KDE::KDE2d kde(NUM1, NUM2,
		     I1min, I1max, I2min, I2max, Sig1, Sig2);

      histoM  = kde(histoM);
      histoC  = kde(histoC);
      histoE  = kde(histoE);
      histoJ  = kde(histoJ);
      histoR  = kde(histoR);
      histoI  = kde(histoI);
      histoT  = kde(histoT);
      histoF  = kde(histoF);
      histo1  = kde(histo1);
      histo2  = kde(histo2);
      histoDF = kde(histoDF);
    }

    for (int j=0; j<NUM2; j++) {

      if (Kpow)
	I2 = pow(pow(I2min, KPOWER) + d3*(0.5+j), 1.0/KPOWER);
      else
	I2 = I2min + d2*(0.5+j);
      for (int i=0; i<NUM1; i++) {
	I1 = I1min + d1*(0.5+i);

	if (not actions) orb.new_orbit(I1, I2);
	if (not actions and relJ) histoJ(i, j) /= orb.Jmax();

	if (histoC(i, j)>MINBIN && histoM(i, j)>0.0) {
	  double jfac = 1.0;
	  if (not actions and jaco)
	    jfac = orb.Jmax()*orb.Jmax()*I2/orb.get_freq(1);
	  if (SPECIFIC) {
	    p_rec(out[0], I1, I2, histoM(i, j)/jfac);
	    p_rec(out[1], I1, I2, histoE(i, j)/histoM(i, j));
	    p_rec(out[2], I1, I2, histoJ(i, j)/histoM(i, j));
	    p_rec(out[3], I1, I2, histoT(i, j)/histoM(i, j));
	    p_rec(out[4], I1, I2, histoR(i, j)/histoM(i, j));
	    p_rec(out[5], I1, I2, histoI(i, j)/histoM(i, j));
	    if (I1avg[i]>0.0)
	      p_rec(out[6], I1, I2, histoJ(i, j)/histoM(i, j)/I1avg[i]);
	    else
	      p_rec(out[6], I1, I2, 0);
	  } else {
	    p_rec(out[0], I1, I2, histoM(i, j)/mfac);
	    p_rec(out[1], I1, I2, histoE(i, j)/mfac);
	    p_rec(out[2], I1, I2, histoJ(i, j)/mfac);
	    p_rec(out[3], I1, I2, histoT(i, j)/mfac);
	    p_rec(out[4], I1, I2, histoR(i, j)/histoM(i, j));
	    p_rec(out[6], I1, I2, histoI(i, j)/histoM(i, j));
	    if (I1avg[i]>0.0)
	      p_rec(out[6], I1, I2, histoJ(i, j)/I1avg[i]/mfac);
	    else
	      p_rec(out[6], I1, I2, 0.0);
	  }
	  p_rec(out[7],  I1, I2, histoM(i, j)/histoC(i, j));
	  p_rec(out[9],  I1, I2, rapo  (i, j)/histoM(i, j));
	  p_rec(out[10], I1, I2, rperi (i, j)/histoM(i, j));
	  p_rec(out[11], I1, I2, omega1(i, j)/histoM(i, j));
	  p_rec(out[12], I1, I2, omega2(i, j)/histoM(i, j));
	} else {
	  for (int k=0; k<8;  k++) p_rec(out[k], I1, I2, 0.0);
	  for (int k=9; k<13; k++) p_rec(out[k], I1, I2, 0.0);
	}

	double jfac = 1.0;
	if (jaco) {
	  if (actions) jfac = 1.0/I2;
	  else         jfac = orb.get_freq(1)/(orb.Jmax()*orb.Jmax()*I2);
	}
	  
	if (totMass>0.0) {
	  p_rec(out[8 ], I1, I2, histoF(i, j)*jfac/totMass);
	  p_rec(out[13], I1, I2, histo1(i, j)*jfac/totMass);
	  p_rec(out[14], I1, I2, histo2(i, j)*jfac/totMass);
	}
	else {
	  p_rec(out[8 ], I1, I2, 0.0);
	  p_rec(out[13], I1, I2, 0.0);
	  p_rec(out[14], I1, I2, 0.0);
	}

	p_rec(out[15], I1, I2, histoDF(i, j));
	p_rec(out[16], I1, I2, histoC (i, j));
      }
      if (not meshgrid) for (int k=0; k<17; k++) out[k] << endl;
    }
    
    if (CUMULATE) {
      string CUMFILE = OUTFILE + ".cum";
      ofstream outc(CUMFILE.c_str(), ios_base::out | ios_base::app);
      if (outc) {
	Eigen::VectorXd delL = Eigen::VectorXd::Zero(NUM1);
	Eigen::VectorXd cumL = I1avg;
	Eigen::VectorXd cumM = I1mass;
	for (int i=0; i<NUM1; i++) {
	  cumL[i] *= I1mass[i];
	  for (int j=0; j<NUM2; j++) delL[i] += histoJ(i, j);
	}
	for (int i=1; i<NUM1; i++) {
	  cumL[i] += cumL[i-1];
	  cumM[i] += cumM[i-1];
	}
	
	for (int i=0; i<NUM1; i++)
	  outc << setw(18) << final_time
	       << setw(18) << I1min+(0.5+i)*d1
	       << setw(18) << delL[i]
	       << setw(18) << ( cumM[i]>0.0 ? delL[i]/cumM[i] : 0.0 )
	       << setw(18) << ( cumL[i]>0.0 ? delL[i]/cumL[i] : 0.0 )
	       << endl;
	
	outc << endl;
      }
    }
    
    const int nrlabs = 7, fieldsz=18;
    const char *rlabels[] = {
			     "Radius",
			     "Delta E",
			     "Delta Lz",
			     "(Delta E)/<E>",
			     "(Delta L)/<L>",
			     "<M>",
			     "<J>"
    };
    
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[17] << "#";
      else      out[17] << "+";
      out[17] << setw(fieldsz-1) << left << setfill('-') << '-';
    }
    out[17] << endl << setfill(' ');
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[17] << "# ";
      else      out[17] << "+ ";
      out[17] << setw(fieldsz-2) << left << rlabels[j];
    }
    out[17] << endl;
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[17] << "# ";
      else      out[17] << "+ ";
      out[17] << setw(fieldsz-2) << left << j+1;
    }
    out[17] << endl;
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[17] << "#";
      else      out[17] << "+";
      out[17] << setw(fieldsz-1) << left << setfill('-') << '-';
    }
    out[17] << endl << setfill(' ');
    
    for (int i=0; i<NUMR; i++) {
      
      double rr = rhmin + dR*(0.5+i);
      if (LOGR) rr = exp(rr);
      
      out[18] << setw(fieldsz) << rr 
	      << setw(fieldsz) << histoP[i]
	      << setw(fieldsz) << histoL[i]
	      << setw(fieldsz) << histPr[i]
	      << setw(fieldsz) << histLr[i]
	      << setw(fieldsz) << histoN[i]
	      << setw(fieldsz) << histoS[i]
	      << endl;
    }
    
    {
      const int nrlabs = 5, fieldsz=18;

      std::vector<std::vector<std::string>>
	labels = { {"E", "F1", "F2", "DF", "Df"}, 
		   {"K", "F1", "F2", "DF", "Df"} };

      if (actions) {
	labels[0][0] = "I_1";
	labels[1][0] = "I_2";
      }

      for (int l=0; l<2; l++) {

	for (int j=0; j<nrlabs; j++) {
	  if (j==0) out[19+l] << "#";
	  else      out[19+l] << "+";
	  out[19+l] << setw(fieldsz-1) << left << setfill('-') << '-';
	}

	out[19+l] << endl << setfill(' ');
	
	for (int j=0; j<nrlabs; j++) {
	  if (j==0) out[19+l] << "# ";
	  else      out[19+l] << "+ ";
	  out[19+l] << setw(fieldsz-2) << left << labels[l][j];
	}
	
	out[19+l] << endl;
	for (int j=0; j<nrlabs; j++) {
	  if (j==0) out[19+l] << "# ";
	  else      out[19+l] << "+ ";
	  out[19+l] << setw(fieldsz-2) << left << j+1;
	}
	out[19+l] << endl;

	for (int j=0; j<nrlabs; j++) {
	  if (j==0) out[19+l] << "#";
	  else      out[19+l] << "+";
	  out[19+l] << setw(fieldsz-1) << left << setfill('-') << '-';
	}
	out[19+l] << endl << setfill(' ');
	
	for (int j=0; j<histo1_1d[l].size(); j++) {
	  if (l==0) out[19+l] << setw(fieldsz) << left << I1min + d1*j;
	  else      out[19+l] << setw(fieldsz) << left << I2min + d2*j;
	  
	  out[19+l] << setw(fieldsz) << left << histo1_1d[l][j]
		    << setw(fieldsz) << left << histo2_1d[l][j]
		    << setw(fieldsz) << left << histoF_1d[l][j]
		    << setw(fieldsz) << left << histoDF_1d[l][j]
		    << std::endl;
	}
      }
    }
  }
  
  MPI_Finalize();
  
  return 0;
}
