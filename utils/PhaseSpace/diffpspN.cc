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

#include <ParticleReader.H>
#include <massmodel.H>
#include <localmpi.H>
#include <EXPini.H>		// Enhanced option parsing
#include <libvars.H>		// EXP library globals
#include <interp.H>

void p_rec(std::ofstream& out, double E, double K, double V)
{
  out << std::setw(15) << E
      << std::setw(15) << K
      << std::setw(15) << V << std::endl;
}


double p1[3] = {0.0, 0.0, 0.0};	// Phase space #1 center
double p2[3] = {0.0, 0.0, 0.0};	// Phase space #2 center

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
  double       MINBIN;
  bool         SPECIFIC;
  bool         CUMULATE;
  bool         LZDIST;
  bool         LOGMOD;
  bool         meshgrid;
  double       DIVERGE_RFAC;
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
    "   OUTFILE.DM       Mass per bin                   [2d]\n"		\
    "   OUTFILE.DE       Delta E                        [2d]\n"		\
    "   OUTFILE.DK       Delta L (tangen action)        [2d]\n"		\
    "   OUTFILE.DI       Delta I (radial action)        [2d]\n"		\
    "   OUTFILE.DJ       Ang mom per bin                [2d]\n"		\
    "   OUTFILE.DKJ      Delta J/J_avg(E)               [2d]\n"		\
    "   OUTFILE.Dm       Mean mass per bin              [2d]\n"		\
    "   OUTFILE.DF       Distribution function          [2d]\n"		\
    "   OUTFILE.DR       Run mass, J, Delta J (R)       [1d]\n"		\
    "   OUTFILE.chk      Orbital element check              \n"		\
    "=======================================================\n"
    " E=Energy, K=J/J_{max}(E)                              \n"
    "=======================================================\n";



  cxxopts::Options options(argv[0], desc);

  options.add_options()
    ("h,help", "This help message")
    ("rmaxf", "Do not extend orbit evaluation beyond the model radius")
    ("relJ", "Compute the relative binned angular momentum")
    ("jaco", "Compute phase-space Jacobian for DF computation")
    ("actions", "Print output in action space rather than E-kappa space.  The default is Energy-Kappa.")
    ("F,filetype", "input file type (one of: PSPout, PSPspl, GadgetNative, GadgetHDF5)",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("RMIN", "Minimum model radius",
     cxxopts::value<double>(RMIN)->default_value("0.0"))
    ("RMAX", "Maximum model radius",
     cxxopts::value<double>(RMAX)->default_value("1.0"))
    ("BMIN", "Minimum value in cos(b)",
     cxxopts::value<double>(BMIN)->default_value("-1.0"))
    ("BMAX", "Maximum value in cos(b)",
     cxxopts::value<double>(BMAX)->default_value("1.0"))
    ("EMIN", "Minimum value in phase-space energy",
      cxxopts::value<double>(EMIN)->default_value("-1.0e20"))
    ("EMAX", "Maximum value in phase-space energy",
     cxxopts::value<double>(EMAX)->default_value("1.0e20"))
    ("KMIN", "Minimum value in J/J_max}(E)",
     cxxopts::value<double>(KMIN)->default_value("0.0"))
    ("KMAX", "Minimum value in J/J_max}(E)",
      cxxopts::value<double>(KMAX)->default_value("1.0"))
    ("NLIMIT", "Number of particles to compare (NLIMIT<0 means all)",
     cxxopts::value<int>(NLIMIT)->default_value("-1"))
    ("NSKIP", "Number of particles to skip before comparing",
     cxxopts::value<int>(NSKIP)->default_value("0"))
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
    ("DIVERGE_RFAC", "Cusp index for extrapolation",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("1.0"))
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
  for (auto file : INFILE1) {
    if (not std::filesystem::exists(std::filesystem::path(file))) bad++;
  }

  if (bad) {
    if (myid==0)
      std::cerr << "Could not open " << bad
		<< " file(s) in INFILE1 list" << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  for (auto file : INFILE2) {
    if (not std::filesystem::exists(std::filesystem::path(file))) bad++;
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
  Eigen::MatrixXd histoC, histoM, histoE, histoJ, histoI, histoT;
  Eigen::MatrixXd histo1, histo2;

  histoC = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoM = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoE = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoJ = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoI = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histoT = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histo1 = Eigen::MatrixXd::Zero(NUM1, NUM2);
  histo2 = Eigen::MatrixXd::Zero(NUM1, NUM2);

  Eigen::VectorXd histoP, histoL, histPr, histLr, histoS, histoN;

  histoP = Eigen::VectorXd::Zero(NUMR);
  histoL = Eigen::VectorXd::Zero(NUMR);
  histPr = Eigen::VectorXd::Zero(NUMR);
  histLr = Eigen::VectorXd::Zero(NUMR);
  histoS = Eigen::VectorXd::Zero(NUMR);
  histoN = Eigen::VectorXd::Zero(NUMR);

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
    
  const int nfiles = 11;
  const char *suffix[] = {
    ".DM", 			// Mass per bin (E, K)		#0
    ".DE", 			// Delta E(E, K)		#1
    ".DK", 			// Delta J(E, K)		#2
    ".DJ", 			// Ang mom per bin (E, K)	#3
    ".DI", 			// Rad. act. per bin (E, K)	#4
    ".DKJ", 			// Delta J(E, K)/J_avg(E)	#5
    ".Dm", 			// Mean mass per bin (E, K)	#6
    ".DF", 			// Delta DF (E, K)              #7
    ".Df", 			// Delta DF/F (E, K)            #8
    ".DR", 			// Run mass, J, Delta J (R)	#9
    ".chk"			// Orbital element check	#10
  };
  std::vector<string> filename(nfiles);
  for (int i=0; i<nfiles; i++) filename[i] = OUTFILE + suffix[i];

  std::vector<std::ofstream> out(nfiles);
  
  bool good = true;
  if (myid==0) {
    for (int i=0; i<nfiles; i++) {
    
      out[i].open(filename[i].c_str());

      if (!out[i]) {
	cerr << "Error opening <" << filename[i] << ">" << std::endl;
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

  PR::PSPstanza *stanza, *stanza1, *stanza2;

  PR::PRptr psp;

  auto hmodel = std::make_shared<SphericalModelTable>(MODELFILE, DIVERGE, DIVERGE_RFAC);

  if (vm.count("rmaxf")) {
    SphericalOrbit::RMAXF = 1.0;
    std::cout << "RMAXF=" << SphericalOrbit::RMAXF << std::endl;
  }

  SphericalOrbit orb(hmodel);

  double Emin = std::max<double>(hmodel->get_pot(hmodel->get_min_radius()), EMIN);
  double Emax = std::min<double>(hmodel->get_pot(hmodel->get_max_radius()), EMAX);

  double I1min, I1max, I2min, I2max;

  if (actions) {
    Emin *= 1.0 - KTOL;
    Emax *= 1.0 + KTOL;

    orb.new_orbit(Emin, 1.0-KTOL);
    I1min = orb.get_action(0);

    orb.new_orbit(Emin, KTOL);
    I2min = orb.get_action(1);

    orb.new_orbit(Emax, KTOL);
    I1max = orb.get_action(0);

    orb.new_orbit(Emax, 1.0 - KTOL);
    I2max = orb.get_action(1);

  } else {
    I1min = Emin;
    I1max = Emax;
    I2min = KTOL;
    I2max = 1.0 - KTOL;
  }

  double d1 = (I1max - I1min) / NUM1;
  double d2 = (I2max - I2min) / NUM2;
  double I1, I2;


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
  int numK0=0, numK1=0, reject=0, N=0;
  double KOVER  = 0.0;
  double KUNDER = 1.0;

  // Times for the PSP snaps
  //
  double initl_time, final_time;


  // Iterate through file list
  //
  for (size_t n=0; n<INFILE1.size(); n++) {

    PR::PRptr psp1, psp2;

    try {
      psp1 = PR::ParticleReader::createReader(fileType, {INFILE1[n]}, myid, true);
  
      initl_time = psp1->CurrentTime();

      psp1->SelectType(COMP);

      if (myid==0) {
	std::cout << std::endl << std::string(40, '-') << std::endl;
	std::cout << "File 1: " << INFILE1[n] << std::endl;
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
      psp2 = PR::ParticleReader::createReader(fileType, {INFILE2[n]}, myid, true);
  
      final_time = psp2->CurrentTime();

      psp2->SelectType(COMP);

      if (myid==0) {
	std::cout << "File 2: " << INFILE2[n] << endl;
      std::cout << "Found dump at time: " << final_time << std::endl;
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
    // Do difference
    //============================================================
    
    struct SPS
    {
      double pos[3], vel[3];
    } tps;
    

    std::map<int, SPS> ph;

    N = 0;
    for (auto pp=psp1->firstParticle(); pp!=0; pp=psp1->nextParticle()) {
      
      if (N++ % numprocs == myid) {
	for (int k=0; k<3; k++) {
	  tps.pos[k] = pp->pos[k];
	  tps.vel[k] = pp->vel[k];
	}
	ph[pp->indx] = tps;
      }
    }

    N = 0;
    for (auto pp=psp2->firstParticle(); pp!=0; pp=psp2->nextParticle()) {
    
      auto ip = ph.find(pp->indx);
    
      if (ip != ph.end()) {
      
	double angmom1[3], angmom2[3];
	double p10[3], p20[3], v10[3], v20[3];
      
	if (TAG>=0 && stanza1->comp.niatr>TAG) {
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
	  E1 = max<double>(E1, Emin);
	  E1 = min<double>(E1, Emax);
	  
	  double E2 = 0.5*vv2 + hmodel->get_pot(sqrt(rr2));
	  E2 = max<double>(E2, Emin);
	  E2 = min<double>(E2, Emax);
	  
	  angmom1[0] = p10[1]*v10[2] - p10[2]*v10[1];
	  angmom1[1] = p10[2]*v10[0] - p10[0]*v10[2];
	  angmom1[2] = p10[0]*v10[1] - p10[1]*v10[0];
	
	  double cosb = angmom1[2]/sqrt(angmom1[0]*angmom1[0] +
					angmom1[1]*angmom1[1] +
					angmom1[2]*angmom1[2] );
	
	  if (POSNEG>0 && angmom1[2]<0.0 || POSNEG<0 && angmom1[2]>0.0)
	    continue;

	  if (cosb<BMIN || cosb>BMAX)
	    continue;
	
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
	    if (K1>1.0-KTOL) {
	      numK1++;
	      KOVER = std::max<double>(KOVER,  K1);
	      K1 = 1.0 - KTOL;
	    }
	    if (K1<KTOL) {
	      numK0++;
	      KUNDER = min<double>(KUNDER, K1);
	      K1 = KTOL;
	    }
	    
	    orb.new_orbit(E1, K1);
	  
	    I1 = orb.get_action(0);
	  
	    orb.new_orbit(E2, 0.5);
	    double K2 = sqrt(j2)/orb.Jmax();
	    if (K2>1.0-KTOL) {
	      numK1++;
	      KOVER  = max<double>(KOVER, K2);
	      K2 = 1.0 - KTOL;
	    }
	    if (K2<KTOL) {
	      numK0++;
	      KUNDER = min<double>(KUNDER, K2);
	      K2 = KTOL;
	    }
	  
	    orb.new_orbit(E2, K2);

	    double I2 = orb.get_action(1);
	  
	    if (myid==0) {
	      if (WHICHEK & 1) {
		if (K1>1.0 || K1<0.0)
		  out[10] << setw(15) << E1 << setw(15) << K1
			  << setw(15) << E2 << setw(15) << sqrt(j1)
			  << setw(15) << sqrt(rr1) << setw(15) << sqrt(rr2)
			  << setw(5) << 1 << endl;
	      }
	    
	      if (WHICHEK & 2) {
		if (K2>1.0 || K2<0.0)
		  out[10] << setw(15) << E2 << setw(15) << K2
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
	    
	    KK = std::min<double>(KK, 1.0-KTOL);
	    KK = std::max<double>(KK, KTOL);
	  
	    double I1_1, I2_1, I1_2, I2_2;
	    int i1, i2, i11, i12, i21, i22;

	    if (actions) {

	      // Get the actions
	      orb.new_orbit(EE, KK);

	      I1 = orb.get_action(0);
	      I2 = orb.get_action(1);

	      I1 = std::max<double>(I1, I1min);
	      I1 = std::min<double>(I1, I1max);
	  
	      I2 = std::max<double>(I2, I2min);
	      I2 = std::min<double>(I2, I2max);
	      
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

	      i11 = (int)floor( (E1 - Emin) / d1 );
	      i11 = std::max<int>(i11, 0);
	      i11 = std::min<int>(i11, NUM1-1);
	    
	      i21 = (int)floor( K1 / d2 );
	      i21 = std::max<int>(i21, 0);
	      i21 = std::min<int>(i21, NUM2-1);

	      i12 = (int)floor( (E2 - Emin) / d1 );
	      i12 = std::max<int>(i12, 0);
	      i12 = std::min<int>(i12, NUM1-1);
	  
	      i22 = (int)floor( K2 / d2 );
	      i22 = max<int>(i22, 0);
	      i22 = min<int>(i22, NUM2-1);

	      i1 = (int)floor( (EE - Emin) / d1 );
	      i1 = max<int>(i1, 0);
	      i1 = min<int>(i1, NUM1-1);
	    
	      i2 = (int)floor( KK / d2 );
	      i2 = max<int>(i2, 0);
	      i2 = min<int>(i2, NUM2-1);
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
	    histoI(i1, i2) += pp->mass*(I2 - I1);
	    histoT(i1, i2) += pp->mass*L1;
	    
	    histo1(i11, i21) += pp->mass;
	    histo2(i12, i22) += pp->mass;
	    
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
	  }
	  catch (const std::runtime_error& error) {
	    reject++;
	  }
	}
      }
      
      if (myid==0 and NREPORT) {
	if (!((N+1)%NREPORT)) cout << "\rProcessed: " 
				   << setw(10) << (N+1)*numprocs << flush;
	N++;
      }
    }

    if (myid==0 and NREPORT)
      std::cout << std::endl << std::string(40, '-') << std::endl;
  }
  
  
  // Send to root
  //
  if (myid) {
    MPI_Reduce(&reject,       0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoC.data(), 0, histoC.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoM.data(), 0, histoM.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo1.data(), 0, histo1.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histo2.data(), 0, histo2.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoE.data(), 0, histoE.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoJ.data(), 0, histoJ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoI.data(), 0, histoI.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoT.data(), 0, histoT.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoP.data(), 0, histoP.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoL.data(), 0, histoL.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histPr.data(), 0, histPr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histLr.data(), 0, histLr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoS.data(), 0, histoS.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoN.data(), 0, histoN.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  // Receive at root
  //
  else {
    MPI_Reduce(MPI_IN_PLACE, &reject, 1,                   MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoC.data(), histoC.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoM.data(), histoM.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo1.data(), histo1.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histo2.data(), histo2.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoE.data(), histoE.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoJ.data(), histoJ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoI.data(), histoI.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoT.data(), histoT.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoP.data(), histoP.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoL.data(), histoL.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histPr.data(), histPr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histLr.data(), histLr.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoS.data(), histoS.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, histoN.data(), histoN.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    std::cout << std::endl;
    
    if (reject) std::cout << "SphericalOrbit failures in " << reject
			  << "/" << N << " states" << std::endl << std::endl;
    
    Eigen::VectorXd I1avg  = Eigen::VectorXd::Zero(NUM1);
    Eigen::VectorXd I2avg  = Eigen::VectorXd::Zero(NUM2);
    Eigen::VectorXd I1mass = Eigen::VectorXd::Zero(NUM1);
    Eigen::VectorXd I2mass = Eigen::VectorXd::Zero(NUM2);

    // Delta DF
    Eigen::MatrixXd histoF  = histo2 - histo1;
    Eigen::MatrixXd histoDF = histoF;
    
    // Relative Delta DF
    for (size_t j=0; j<histoF.cols(); j++) { // Footron order...
      for (size_t i=0; i<histoF.rows(); i++) {
	if (histo1(i, j)>0.0) histoDF(i, j) = histoF(i, j)/histo1(i, j);
	else                  histoDF(i, j) = 0.0;
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
      for (int k=0; k<9; k++) out[k] << std::setw(8) << NUM1
				     << std::setw(8) << NUM2
				     << std::endl;
    bool relJ = false;
    if (vm.count("relJ")) relJ = true;

    for (int j=0; j<NUM2; j++) {
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
	    p_rec(out[4], I1, I2, histoI(i, j)/histoM(i, j));
	    if (I1avg[i]>0.0)
	      p_rec(out[5], I1, I2, histoJ(i, j)/histoM(i, j)/I1avg[i]);
	    else
	      p_rec(out[5], I1, I2, 0);
	  } else {
	    p_rec(out[0], I1, I2, histoM(i, j)/mfac);
	    p_rec(out[1], I1, I2, histoE(i, j)/mfac);
	    p_rec(out[2], I1, I2, histoJ(i, j)/mfac);
	    p_rec(out[3], I1, I2, histoT(i, j)/mfac);
	    p_rec(out[4], I1, I2, histoI(i, j)/histoM(i, j));
	    if (I1avg[i]>0.0)
	      p_rec(out[5], I1, I2, histoJ(i, j)/I1avg[i]/mfac);
	    else
	      p_rec(out[5], I1, I2, 0.0);
	  }
	  p_rec(out[6], I1, I2, histoM(i, j)/histoC(i, j));
	} else {
	  for (int k=0; k<7; k++) p_rec(out[k], I1, I2, 0.0);
	}

	double jfac = 1.0;
	if (jaco) {
	  if (actions) jfac = 1.0/I2;
	  else         jfac = orb.get_freq(1)/(orb.Jmax()*orb.Jmax()*I2);
	}
	  
	if (totMass>0.0) {
	  p_rec(out[7], I1, I2, histoF(i, j)*jfac/totMass);
	}
	else {
	  p_rec(out[7], I1, I2, 0.0);
	}

	p_rec(out[8], I1, I2, histoDF(i, j));
      }
      if (not meshgrid) for (int k=0; k<9; k++) out[k] << endl;
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
      if (j==0) out[9] << "#";
      else      out[9] << "+";
      out[9] << setw(fieldsz-1) << left << setfill('-') << '-';
    }
    out[9] << endl << setfill(' ');
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[9] << "# ";
      else      out[9] << "+ ";
      out[9] << setw(fieldsz-2) << left << rlabels[j];
    }
    out[9] << endl;
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[9] << "# ";
      else      out[9] << "+ ";
      out[9] << setw(fieldsz-2) << left << j+1;
    }
    out[9] << endl;
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[9] << "#";
      else      out[9] << "+";
      out[9] << setw(fieldsz-1) << left << setfill('-') << '-';
    }
    out[9] << endl << setfill(' ');
    
    for (int i=0; i<NUMR; i++) {
      
      double rr = rhmin + dR*(0.5+i);
      if (LOGR) rr = exp(rr);
      
      out[9] << setw(fieldsz) << rr 
	     << setw(fieldsz) << histoP[i]
	     << setw(fieldsz) << histoL[i]
	     << setw(fieldsz) << histPr[i]
	     << setw(fieldsz) << histLr[i]
	     << setw(fieldsz) << histoN[i]
	     << setw(fieldsz) << histoS[i]
	     << endl;
    }
    
    cout << "K check:" << endl
	 << "  Under=" << numK0 << "  Min=" << KUNDER << endl
	 << "   Over=" << numK1 << "  Max=" << KOVER  << endl;
  }
  
  MPI_Finalize();
  
  return 0;
}
