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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <Eigen/Eigen>

#include <ParticleReader.H>
#include <massmodel.H>
#include <localmpi.H>
#include <cxxopts.H>		// Option parsing
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
  const double KTOL = 1.0e-2;

  // MPI initialization
  //
  local_init_mpi(argc, argv);

  // Parameter assignment
  //
  double       TIME1;
  double       TIME2;
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
  int          NUME;
  int          NUMK;
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

  std::vector<std::string> INFILE1, INFILE2;

  const char* desc = 
    "=======================================================\n"		\
    "Compute resonance loci for a given peturbing frequency \n"		\
    "=======================================================\n"		\
    "   Output file key:\n"						\
    "   ----------------\n"						\
    "   OUTFILE.DM       Mass per bin (E, K)            [2d]\n"		\
    "   OUTFILE.DE       Delta E(E, K)                  [2d]\n"		\
    "   OUTFILE.DK       Delta H(E, K) (tangen action)  [2d]\n"		\
    "   OUTFILE.DI       Delta I(E, K) (radial action)  [2d]\n"		\
    "   OUTFILE.DJ       Ang mom per bin (E, K)         [2d]\n"		\
    "   OUTFILE.DKJ      Delta J(E, K)/J_avg(E)         [2d]\n"		\
    "   OUTFILE.Dm       Mean mass per bin (E, K)       [2d]\n"		\
    "   OUTFILE.DF       Distribution function (E, K)   [2d]\n"		\
    "   OUTFILE.DR       Run mass, J, Delta J (R)       [1d]\n"		\
    "   OUTFILE.chk      Orbital element check              \n"		\
    "=======================================================\n"
    " E=Energy, K=J/J_{max}(E)                              \n"
    "=======================================================\n";



  cxxopts::Options options(argv[0], desc);

  options.add_options()
    ("h,help", "This help message")
    ("2,kappa2", "Bin kappa^2 instead of kappa")
    ("rmaxf", "Do not extend orbit evaluation beyond the model radius")
    ("relJ", "Compute the relative binned angular momentum")
    ("jaco", "Compute phase-space Jacobian for DF computation")
    ("F,filetype", "input file type (one of: PSPout, PSPspl, GadgetNative, GadgetHDF5)",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("TIME1", "Target time for fiducial phase space",
     cxxopts::value<double>(TIME1)->default_value("0.0"))
    ("TIME2", "Target time for evolved phase space",
     cxxopts::value<double>(TIME2)->default_value("10.0"))
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
    ("NUME", "Number of energy bins",
     cxxopts::value<int>(NUME)->default_value("60"))
    ("NUMK", "Number of kappa bins",
     cxxopts::value<int>(NUMK)->default_value("20"))
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
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0) std::cout << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }

  if (INFILE1.size() != INFILE2.size()) {
    if (myid==0)
      std::cerr << "Input file vectors must have paired entries"
		<< std::endl;
    MPI_Finalize();
    exit(-1);
  }

  bool kappa2 = false;
  if (vm.count("kappa2")) kappa2 = true;

  std::ofstream tout;
  if (myid==0) tout.open("tmp.ktest");

  bool jaco = false;
  if (vm.count("jaco")) jaco = true;

  // Allocate histograms
  //
  Eigen::MatrixXd histoC, histoM, histoE, histoJ, histoI, histoT;
  Eigen::MatrixXd histo1, histo2;

  histoC = Eigen::MatrixXd::Zero(NUME, NUMK);
  histoM = Eigen::MatrixXd::Zero(NUME, NUMK);
  histoE = Eigen::MatrixXd::Zero(NUME, NUMK);
  histoJ = Eigen::MatrixXd::Zero(NUME, NUMK);
  histoI = Eigen::MatrixXd::Zero(NUME, NUMK);
  histoT = Eigen::MatrixXd::Zero(NUME, NUMK);
  histo1 = Eigen::MatrixXd::Zero(NUME, NUMK);
  histo2 = Eigen::MatrixXd::Zero(NUME, NUMK);

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

  double Emin = max<double>(hmodel->get_pot(hmodel->get_min_radius()), EMIN);
  double Emax = min<double>(hmodel->get_pot(hmodel->get_max_radius()), EMAX);
  double delE = Emax - Emin;
  double dK = (1.0 - KTOL)  / NUMK;
  double dE = delE / NUME;

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
      psp1 = PR::ParticleReader::createReader(fileType, INFILE1[n], myid, true);
  
      initl_time = psp1->CurrentTime();

      psp1->SelectType(COMP);

      if (myid==0) {
	std::cout << "File 1: " << INFILE1[0] << std::endl;
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
      psp2 = PR::ParticleReader::createReader(fileType, INFILE2[n], myid, true);
  
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
	    if (kappa2) K1 *= K1;
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
	    
	    if (kappa2)
	      orb.new_orbit(E1, sqrt(K1));
	    else
	      orb.new_orbit(E1, K1);
	  
	    double I1 = orb.get_action(1);
	  
	    orb.new_orbit(E2, 0.5);
	    double K2 = sqrt(j2)/orb.Jmax();
	    if (kappa2) K2 *= K2;
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
	  
	    if (kappa2)
	      orb.new_orbit(E2, sqrt(K2));
	    else
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
	    
	    tout << KK << endl;
	    
	    KK = min<double>(KK, 1.0);
	    KK = max<double>(KK, 0.0);
	  
	    int ie1 = (int)floor( (E1 - Emin) / dE );
	    ie1 = max<int>(ie1, 0);
	    ie1 = min<int>(ie1, NUME-1);
	    
	    int ik1 = (int)floor( K1 / dK );
	    ik1 = max<int>(ik1, 0);
	    ik1 = min<int>(ik1, NUMK-1);

	    int ie2 = (int)floor( (E2 - Emin) / dE );
	    ie2 = max<int>(ie2, 0);
	    ie2 = min<int>(ie2, NUME-1);
	  
	    int ik2 = (int)floor( K2 / dK );
	    ik2 = max<int>(ik2, 0);
	    ik2 = min<int>(ik2, NUMK-1);

	    int ie = (int)floor( (EE - Emin) / dE );
	    ie = max<int>(ie, 0);
	    ie = min<int>(ie, NUME-1);
	    
	    int ik = (int)floor( KK / dK );
	    ik = max<int>(ik, 0);
	    ik = min<int>(ik, NUMK-1);
	    
	    double L1 = 0.0, L2 = 0.0;
	    for (int k=0; k<3; k++) {
	      L1 += angmom1[k]*angmom1[k];
	      L2 += angmom2[k]*angmom2[k];
	    }
	    L1 = sqrt(L1);
	    L2 = sqrt(L2);
	  
	    histoC(ie, ik) += 1;
	    histoM(ie, ik) += pp->mass;
	    histoE(ie, ik) += pp->mass*(E1 - E2);
	    histoJ(ie, ik) += pp->mass*(L1 - L2);
	    histoI(ie, ik) += pp->mass*(I2 - I1);
	    histoT(ie, ik) += pp->mass*L1;
	    
	    histo1(ie1, ik1) += pp->mass;
	    histo2(ie2, ik2) += pp->mass;
	    
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
				   << setw(10) << N+1 << flush;
	N++;
      }
    }
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
    
    Eigen::VectorXd Eavg  = Eigen::VectorXd::Zero(NUME);
    Eigen::VectorXd Emass = Eigen::VectorXd::Zero(NUME);

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
    
    for (int i=0; i<NUME; i++) {
      
      for (int j=0; j<NUMK; j++) {
	Eavg[i]  += histoT(i, j);
	Emass[i] += histoM(i, j);
	totMass  += histoM(i, j);
      }
      if (Emass[i] > 0.0) Eavg[i] /= Emass[i];
      else                Eavg[i]  = 0.0;
    }
    
    double mfac = dE*dK;
    
    if (meshgrid)
      for (int k=0; k<9; k++) out[k] << std::setw(8) << NUME
				     << std::setw(8) << NUMK
				     << std::endl;
    bool relJ = false;
    if (vm.count("relJ")) relJ = true;

    for (int j=0; j<NUMK; j++) {
      double KK = KTOL + 0.5*dK + dK*j;
      for (int i=0; i<NUME; i++) {
	double EE = Emin + 0.5*dE + dE*i;

	orb.new_orbit(EE, KK);
	if (relJ) histoJ(i, j) /= orb.Jmax();

	if (histoC(i, j)>MINBIN && histoM(i, j)>0.0) {
	  double jfac = 1.0;
	  if (jaco) jfac = orb.Jmax()*orb.Jmax()*KK/orb.get_freq(1);
	  if (SPECIFIC) {
	    p_rec(out[0], EE, KK, histoM(i, j)/jfac);
	    p_rec(out[1], EE, KK, histoE(i, j)/histoM(i, j));
	    p_rec(out[2], EE, KK, histoJ(i, j)/histoM(i, j));
	    p_rec(out[3], EE, KK, histoT(i, j)/histoM(i, j));
	    p_rec(out[4], EE, KK, histoI(i, j)/histoM(i, j));
	    if (Eavg[i]>0.0)
	      p_rec(out[5], EE, KK, histoJ(i, j)/histoM(i, j)/Eavg[i]);
	    else
	      p_rec(out[5], EE, KK, 0);
	  } else {
	    p_rec(out[0], EE, KK, histoM(i, j)/jfac/mfac);
	    p_rec(out[1], EE, KK, histoE(i, j)/mfac);
	    p_rec(out[2], EE, KK, histoJ(i, j)/mfac);
	    p_rec(out[3], EE, KK, histoT(i, j)/mfac);
	    p_rec(out[4], EE, KK, histoI(i, j)/histoM(i, j));
	    if (Eavg[i]>0.0)
	      p_rec(out[5], EE, KK, histoJ(i, j)/Eavg[i]/mfac);
	    else
	      p_rec(out[5], EE, KK, 0.0);
	  }
	  p_rec(out[6], EE, KK, histoM(i, j)/histoC(i, j));
	} else {
	  for (int k=0; k<7; k++) p_rec(out[k], EE, KK, 0.0);
	}

	if (totMass>0.0) {
	  double jfac = 1.0;
	  if (jaco) jfac = orb.Jmax()*orb.Jmax()*KK/orb.get_freq(1);
	  p_rec(out[7], EE, KK, histoF(i, j)/totMass/jfac);
	}
	else {
	  p_rec(out[7], EE, KK, 0.0);
	}

	p_rec(out[8], EE, KK, histoDF(i, j));
      }
      if (not meshgrid) for (int k=0; k<9; k++) out[k] << endl;
    }
    
    if (CUMULATE) {
      string CUMFILE = OUTFILE + ".cum";
      ofstream outc(CUMFILE.c_str(), ios_base::out | ios_base::app);
      if (outc) {
	Eigen::VectorXd delL = Eigen::VectorXd::Zero(NUME);
	Eigen::VectorXd cumL = Eavg;
	Eigen::VectorXd cumM = Emass;
	for (int i=0; i<NUME; i++) {
	  cumL[i] *= Emass[i];
	  for (int j=0; j<NUMK; j++) delL[i] += histoJ(i, j);
	}
	for (int i=1; i<NUME; i++) {
	  cumL[i] += cumL[i-1];
	  cumM[i] += cumM[i-1];
	}
	
	for (int i=0; i<NUME; i++)
	  outc << setw(18) << final_time
	       << setw(18) << Emin+(0.5+i)*dE
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
