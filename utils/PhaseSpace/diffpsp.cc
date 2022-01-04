/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Use ParticleReader class to compute phase space distribution
 *  difference
 *
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
 *  MDW 11/20/91, updated 10/22/21
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
#include <MakeModel.H>
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
  int          WHICHPS;
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
  std::string  INFILE1;
  std::string  INFILE2;
  std::string  CURDIR;
  std::string  MODELFILE;
  std::string  ORIENTFILE;
  std::string  OUTFILE;
  std::string  fileType;

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
    ("WHICHPS", "Which phase-space file to use (1=first, 2=second)",
     cxxopts::value<int>(WHICHPS)->default_value("0"))
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
     cxxopts::value<string>(INFILE1)->default_value("test.bin"))
    ("INFILE2", "Evolved phase-space file",
     cxxopts::value<string>(INFILE2)->default_value("test.bin"))
    ("CURDIR", "Alternative directory",
     cxxopts::value<string>(CURDIR)->default_value(""))
    ("MODELFILE", "Model file for phase-space orbital values",
     cxxopts::value<string>(MODELFILE)->default_value("SLGridSph.model"))
    ("COMP", "Compute wake for this component name",
     cxxopts::value<std::string>(COMP)->default_value("stars"))
    ("ORIENTFILE", "EXP generated orient file for center selection (ignored if null)",
     cxxopts::value<std::string>(ORIENTFILE)->default_value(""))
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

  bool kappa2 = false;
  if (vm.count("kappa2")) kappa2 = true;

  std::ofstream tout;
  if (myid==0) tout.open("tmp.ktest");

  if (WHICHEK < 1 or WHICHEK > 3) {
    if (myid==0)
      std::cerr << "Invalid value for WHICHEK (" << WHICHEK
		<< "): must be 1, 2, or 3\n";
    MPI_Finalize();
    exit(-1);
  }

  //
  // Allocate histograms
  //

  Eigen::MatrixXd histoC, histoM, histoE, histoJ, histoI, histoT;

  histoC = Eigen::MatrixXd::Zero(NUME, NUME);
  histoM = Eigen::MatrixXd::Zero(NUME, NUME);
  histoE = Eigen::MatrixXd::Zero(NUME, NUME);
  histoJ = Eigen::MatrixXd::Zero(NUME, NUME);
  histoI = Eigen::MatrixXd::Zero(NUME, NUME);
  histoT = Eigen::MatrixXd::Zero(NUME, NUME);

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
    
  const int nfiles = 9;
  const char *suffix[] = {
    ".DM", 			// Mass per bin (E, K)		#0
    ".DE", 			// Delta E(E, K)		#1
    ".DK", 			// Delta J(E, K)		#2
    ".DJ", 			// Ang mom per bin (E, K)	#3
    ".DI", 			// Rad. act. per bin (E, K)	#4
    ".DKJ", 			// Delta J(E, K)/J_avg(E)	#5
    ".Dm", 			// Mean mass per bin (E, K)	#6
    ".DR", 			// Run mass, J, Delta J (R)	#7
    ".chk"			// Orbital element check	#8
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

  // ===================================================================
  // Open PSP files
  // ===================================================================

  if (!std::ifstream(INFILE1)) {
    std::cerr << "Error opening <" << INFILE1 << ">" << std::endl;
    exit(-1);
  }

  PR::PRptr psp1, psp2;
  double initl_time, final_time;

  try {
    psp1 = PR::ParticleReader::createReader(fileType, INFILE1, myid, true);
  
    initl_time = psp1->CurrentTime();

    psp1->SelectType(COMP);

    if (myid==0) {
      std::cout << "File 1: " << INFILE1 << std::endl;
      std::cout << "Found dump at time: " << initl_time << std::endl;
    }
  }
  catch (const std::runtime_error& error) {
    if (myid==0)
      std::cerr << "diffpsp: error opening snapshot in file <" << INFILE1 << ">"
		<< std::endl
		<< "diffpsp: " << error.what()
		<< std::endl;
    MPI_Finalize();
    exit(-1);
  }

  try {
    psp2 = PR::ParticleReader::createReader(fileType, INFILE2, myid, true);
  
    final_time = psp2->CurrentTime();

    psp2->SelectType(COMP);

    if (myid==0) {
      std::cout << "File 2: " << INFILE2 << endl;
      std::cout << "Found dump at time: " << final_time << std::endl;
    }
  }
  catch (const std::runtime_error& error) {
    if (myid==0)
      std::cerr << "diffpsp: error opening snapshot in file <" << INFILE2 << ">"
		<< std::endl
		<< "diffpsp: " << error.what()
		<< std::endl;
    MPI_Finalize();
    exit(-1);
  }
    
  // ===================================================================
  // Use orient file?
  // ===================================================================
  
  bool orient = false;
  vector<double> or_time, or_c[3];

  if (ORIENTFILE.length()) {

    orient = true;

    ifstream in(ORIENTFILE.c_str());
    if (!in) {
      if (myid==0)
	std::cerr << "Couldn't open desired orient file: " << ORIENTFILE
		  << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    double val;
    static int buf_size = 1024;
    char buf[buf_size];

    while (!in.eof()) {
      in.getline(buf, buf_size);
      if (in.eof()) break;

      istringstream ins(buf);
      ins >> val;
      or_time.push_back(val);
      for (int k=0; k<5; k++) ins >> val;
      for (int k=0; k<3; k++) {
	ins >> val;
	or_c[k].push_back(val);
      }
      for (int k=0; k<3; k++) {
	ins >> val;
	*(or_c[k].end()-1) += val;
      }
    }
  }
  

  //============================================================
  // Build model
  //============================================================

  std::shared_ptr<SphericalModelTable> hmodel;
  std::shared_ptr<MakeModel> cmodel;
  PR::PSPstanza *stanza, *stanza1, *stanza2;

  PR::PRptr psp;

  if (WHICHPS==0) {
    // Okay, just checking
  }
  else if (WHICHPS==1) {
    psp = psp1;
    std::cout << "Computing model from File #1" << std::endl;
  }
  else if (WHICHPS==2) {
    psp = psp2;
    std::cout << "Computing model from File #2" << std::endl;
  }
  else {
    std::cerr << "Must specify 0, 1 or 2" << std::endl;
    exit(-1);
  }

  if (WHICHPS==0) {

    hmodel = std::make_shared<SphericalModelTable>(MODELFILE, DIVERGE, DIVERGE_RFAC);

    if (vm.count("rmaxf")) {
      SphericalOrbit::RMAXF = 1.0;
      std::cout << "RMAXF=" << SphericalOrbit::RMAXF << std::endl;
    }
  }
  else {

    cmodel = std::make_shared<MakeModel>(RNUM, RMIN, RMAX, LOGMOD);

    int N = 0;

    for (auto part=psp->firstParticle(); part!=0; part=psp->nextParticle()) {

      if (N++ % numprocs == myid) {

	double r = 0.0; 
	for (int k=0; k<3; k++)
	  r += (part->pos[k] - p1[k])*(part->pos[k] - p1[k]);
	r = sqrt(r);
	
	cmodel->AddPoint(r, part->mass);
      }

      if (myid==0 and NREPORT) {
	if (!((N+1)%NREPORT)) std::cout << "\rProcessed: " 
					<< std::setw(10) << N+1 << std::flush;
      }
    }

    
    hmodel = cmodel->Compute();

    if (myid==0) {
      string outfile = OUTFILE + ".file";
      cmodel->WriteModel(outfile);
    }
  }

  SphericalOrbit orb(hmodel);

  double Emin = max<double>(hmodel->get_pot(hmodel->get_min_radius()), EMIN);
  double Emax = min<double>(hmodel->get_pot(hmodel->get_max_radius()), EMAX);
  double delE = Emax - Emin;
  double dK = 1.0  / (NUMK-1);
  double dE = delE / (NUME-1);

  //=======================
  // Open distribution file
  //=======================

  ofstream dout;
  if (LZDIST and myid==0) {
    string DSTFILE = OUTFILE + ".dist";
    dout.open(DSTFILE.c_str());
    if (!dout) {
      cerr << "Couldn't open <" << DSTFILE << ">" << std::endl;
      LZDIST = false;
    }
  }

  //================================
  // Set center of density/potential
  //================================

  if (orient) {
    
    if (psp1->CurrentTime() < or_time.front())
      for (int k=0; k<3; k++) p1[k] = or_c[k].front();
    else if (psp1->CurrentTime() > or_time.back())
      for (int k=0; k<3; k++) p1[k] = or_c[k].back();
    else
      for (int k=0; k<3; k++) p1[k] = odd2(psp1->CurrentTime(), or_time, or_c[k]);

    if (psp2->CurrentTime() < or_time.front())
      for (int k=0; k<3; k++) p2[k] = or_c[k].front();
    else if (psp2->CurrentTime() > or_time.back())
      for (int k=0; k<3; k++) p2[k] = or_c[k].back();
    else
      for (int k=0; k<3; k++) p2[k] = odd2(psp2->CurrentTime(), or_time, or_c[k]);
    
    if (myid==0) {
      cout << setfill('-') << left << setw(70) << "--- Orient center estimates "
	   << endl << setfill(' ') << right;
      cout << "Center at T=" << setw(10) << fixed << psp1->CurrentTime() << ": ";
      for (int k=0; k<3; k++) cout << setw(10) << p1[k];
      cout << endl;
      
      cout << "Center at T=" << setw(10) << fixed << psp2->CurrentTime() << ": ";
      for (int k=0; k<3; k++) cout << setw(10) << p2[k];
      cout << endl;
      cout << setfill('-') << setw(70) << '-' << endl << setfill(' ');
    }
  }

  //============================================================
  // Do difference
  //============================================================

  int numK0=0, numK1=0, reject=0, N=0;
  const double KTOL = 1.0e-2;
  double KOVER  = 0.0;
  double KUNDER = 1.0;

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
	
	if (POSNEG>0 && angmom1[2]<0.0 || POSNEG<0 && angmom1[2]>0.0) continue;
	if (cosb<BMIN || cosb>BMAX) continue;
	
	angmom2[0] = p20[1]*v20[2] - p20[2]*p20[1];
	angmom2[1] = p20[2]*v20[0] - p20[0]*p20[2];
	angmom2[2] = p20[0]*v20[1] - p20[1]*p20[0];
	
	double jj = 0.0, dj = 0.0;
	for (int k=0; k<3; k++) {
	  if (WHICHEK==1)
	    jj += angmom1[k]*angmom1[k];
	  else if (WHICHEK==2)
	    jj += angmom2[k]*angmom2[k];
	  else
	    jj += 0.5*(angmom1[k]*angmom1[k] + angmom2[k]*angmom2[k]);
	  
	  dj += (angmom1[k] - angmom2[k])*(angmom1[k] - angmom2[k]);
	}
	
	try {
	  orb.new_orbit(E1, 0.5);
	  double K1 = sqrt(jj)/orb.Jmax();
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
	  double K2 = sqrt(jj)/orb.Jmax();
	  if (kappa2) K2 *= K2;
	  if (K2>1.0) {
	    numK1++;
	    KOVER  = max<double>(KOVER, K2);
	    K2 = 1.0 - KTOL;
	  }
	  if (K2<0.0) {
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
		out[8] << setw(15) << E1 << setw(15) << K1
		       << setw(15) << E2 << setw(15) << sqrt(jj)
		       << setw(15) << sqrt(rr1) << setw(15) << sqrt(rr2)
		       << setw(5) << 1 << endl;
	    }
	    
	    if (WHICHEK & 2) {
	      if (K2>1.0 || K2<0.0)
		out[8] << setw(15) << E2 << setw(15) << K2
		       << setw(15) << E1 << setw(15) << sqrt(jj)
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
  
  
  // Send to root
  //
  if (myid) {
    MPI_Reduce(&reject,       0,             1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoC.data(), 0, histoC.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(histoM.data(), 0, histoM.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
    
    for (int i=0; i<NUME; i++) {
      
      for (int j=0; j<NUMK; j++) {
	Eavg[i]  += histoT(i, j);
	Emass[i] += histoM(i, j);
      }
      if (Emass[i] > 0.0) Eavg[i] /= Emass[i];
      else                Eavg[i]  = 0.0;
    }
    
    double mfac = dE*dK;
    
    if (meshgrid)
      for (int k=0; k<6; k++) out[k] << std::setw(8) << NUME
				     << std::setw(8) << NUMK
				     << std::endl;
    for (int j=0; j<NUMK; j++) {
      double KK = 0.5*dK + dK*j;
      for (int i=0; i<NUME; i++) {
	double EE = Emin + 0.5*dE + dE*i;
	if (histoC(i, j)>MINBIN && histoM(i, j)>0.0) {
	  if (SPECIFIC) {
	    p_rec(out[0], EE, KK, histoM(i, j));
	    p_rec(out[1], EE, KK, histoE(i, j)/histoM(i, j));
	    p_rec(out[2], EE, KK, histoJ(i, j)/histoM(i, j));
	    p_rec(out[3], EE, KK, histoT(i, j)/histoM(i, j));
	    p_rec(out[4], EE, KK, histoI(i, j)/histoM(i, j));
	    if (Eavg[i]>0.0)
	      p_rec(out[5], EE, KK, histoJ(i, j)/histoM(i, j)/Eavg[i]);
	    else
	      p_rec(out[5], EE, KK, 0);
	  } else {
	    p_rec(out[0], EE, KK, histoM(i, j)/mfac);
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
      }
      if (not meshgrid) for (int k=0; k<7; k++) out[k] << endl;
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
      if (j==0) out[7] << "#";
      else      out[7] << "+";
      out[7] << setw(fieldsz-1) << left << setfill('-') << '-';
    }
    out[7] << endl << setfill(' ');
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[7] << "# ";
      else      out[7] << "+ ";
      out[7] << setw(fieldsz-2) << left << rlabels[j];
    }
    out[7] << endl;
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[7] << "# ";
      else      out[7] << "+ ";
      out[7] << setw(fieldsz-2) << left << j+1;
    }
    out[7] << endl;
    for (int j=0; j<nrlabs; j++) {
      if (j==0) out[7] << "#";
      else      out[7] << "+";
      out[7] << setw(fieldsz-1) << left << setfill('-') << '-';
    }
    out[7] << endl << setfill(' ');
    
    for (int i=0; i<NUMR; i++) {
      
      double rr = rhmin + dR*(0.5+i);
      if (LOGR) rr = exp(rr);
      
      out[7] << setw(fieldsz) << rr 
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
