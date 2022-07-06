#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <cmath>
#include <chrono>
#include <ctime>

using namespace std;

#include "Timer.H"
#include "global.H"
#include "pHOT.H"
#include "TreeDSMC.H"
#include "Collide.H"

#ifdef USE_GPTL
#include <gptl.h>
#endif

#include <signal.h>


static bool CDEBUG      = false; // Thread diagnostics, false for
				 // production

static bool CTIMER      = false; // Thread timing diagnostics, false
				 // for production

bool Collide::PULLIN    = false; // Use the original Pullin velocity
				 // selection algorithm

bool Collide::SampleCRM = false; // Use sample cell for value of mean
				 // collision velocity

// Default NTC quantile threshold
//
double Collide::ntcThreshDef = 0.95;


// Debugging cell length--MFP ratio
//
bool Collide::MFPCL    = true;

// Use the explicit energy solution
//
bool Collide::ESOL     = false;

// Print out sorted cell parameters
//
bool Collide::SORTED   = false;

// Print out T-rho plane for cells 
// with mass weighting
//
bool Collide::PHASE    = false;

// Extra debugging output
//
bool Collide::EXTRA    = false;

// Turn off collisions for testing
//
bool Collide::DRYRUN   = false;

// Turn off cooling for testing
//
bool Collide::NOCOOL   = false;

// Ensemble-based excess cooling
//
bool Collide::ENSEXES  = true;

// Time step diagnostics
//
bool Collide::TSDIAG   = false;

// Cell-volume diagnostics
//
bool Collide::VOLDIAG  = false;

// Mean free path diagnostics
//
bool Collide::MFPDIAG  = false;

// Sample based on maximum (true) or estimate
// from variance (false);
//
bool Collide::NTC      = false;

// Sample based on maximum (true) or estimate
// from variance (false) with no db
//
bool Collide::NTCnodb  = true;

// Use cpu work to augment per particle effort
//
bool Collide::EFFORT   = true;	

// Verbose timing
//
bool Collide::TIMING   = true;

// Temperature floor in EPSM
//
double Collide::TFLOOR = 1000.0;

// Enable mean-free-path time step estimation
//
bool Collide::MFPTS    = false;

// Random number seed
unsigned Collide::seed = 11;

double 				// Enhance (or suppress) fiducial cooling rate
Collide::ENHANCE       = 1.0;	// Currently, only used in LTE method

// Target number of collisions per cell
//
unsigned Collide::collTnum = std::numeric_limits<unsigned>::max();

// Power of two interval for KE/cool histogram
//
int Collide::TSPOW     = 4;

// Proton mass (g)
//
const double mp        = 1.67262158e-24;

// Boltzmann constant (cgs)
//
const double boltz     = 1.3810e-16;

// Allow counter to stop job
//
bool Collide::numSanityStop     = false;

// Maximum number per step
//
unsigned Collide::numSanityMax  = 100000000u;

// Verbose messaging
//
bool Collide::numSanityMsg      = false;

// Lower thresh for reporting
//
unsigned Collide::numSanityVal  = 10000000u;

// Upper thresh for reporting
//
unsigned Collide::numSanityFreq = 2000000u;

// Default interation type
//
std::string Collide::Labels::def = "default";

// Electron type tag
//
const unsigned short Collide::electronTyp = -2;

// Enable NTC diagnostics
//
bool Collide::DEBUG_NTC         = false;

// Enable NTC full distribution
//
bool Collide::NTC_DIST          = false;

// Oversampling factor for TESTSPREAD
//
int Collide::TestSpreadCount    = 40;

// Set cuda
bool Collide::use_cuda          = false;

// For logging output
//
std::string Collide::printDivider(80, '-');

extern "C"
void *
collide_thread_call(void *atp)
{
  thrd_pass_Collide *tp = (thrd_pass_Collide *)atp;
  Collide *p = (Collide *)tp->p;
  p -> collide_thread((void*)&tp->arg);
  return NULL;
}

void Collide::collide_thread_fork(sKeyDmap* Fn)
{
  std::clock_t startcputime;
  std::chrono::high_resolution_clock::time_point wcts;

  if (CTIMER) {
    startcputime = std::clock();
    wcts = std::chrono::high_resolution_clock::now();
  }

  int errcode;
  void *retval;
  
#if HAVE_LIBCUDA==1
  if (nthrds==1 or c0->cudaDevice>=0) {
#else
  if (nthrds==1) {
#endif
    thrd_pass_Collide td;
    
    td.p        = this;
    td.arg.fn   = Fn;
    td.arg.id   = 0;
    
    collide_thread_call(&td);
    
    return;
  }
  
  td = new thrd_pass_Collide [nthrds];
  if (!td) {
    std::ostringstream sout;
    sout << "Process " << myid 
	 << ": Collide::collide_thread_fork: error allocating memory for thread counters\n";
    throw std::runtime_error(sout.str());
  }

  t = new pthread_t [nthrds];
  if (!t) {
    std::ostringstream sout;
    sout << "Process " << myid
	 << ": collide_thread_fork: error allocating memory for thread\n";
    throw std::runtime_error(sout.str());
  }
  
  // Make the <nthrds> threads
  for (int i=0; i<nthrds; i++) {
    td[i].p        = this;
    td[i].arg.fn   = Fn;
    td[i].arg.id   = i;
    
    errcode =  pthread_create(&t[i], 0, collide_thread_call, &td[i]);
    if (errcode) {
      cerr << "Process " << myid;
      cerr << " collide: cannot make thread " << i
	   << ", errcode=" << errcode << endl;
      exit(19);
    }
  }
  
  if (CDEBUG)
    cerr << "Process " << myid << ": " << nthrds << " threads created"
	 << std::endl;
  
  waitTime.start();
  
  // Collapse the threads
  for (int i=0; i<nthrds; i++) {
    if ((errcode=pthread_join(t[i], &retval))) {
      cerr << "Process " << myid;
      cerr << " collide: thread join " << i
	   << " failed, errcode=" << errcode << endl;
      exit(20);
    }
    if (i==0) {
      waitSoFar = waitTime.stop();
      joinTime.start();
    }
  }
  
  if (CDEBUG)
    cerr << "Process " << myid << ": " << nthrds << " threads joined"
	 << std::endl;
  
  
  joinSoFar = joinTime.stop();
  
  delete [] td;
  delete [] t;

  if (CTIMER) {

    double cpu_duration =
      (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;

    std::chrono::duration<double> wctduration =
      (std::chrono::high_resolution_clock::now() - wcts);

    //  +--- Summary output
    //  |
    //  V
    if (true) {
      if (wctduration > std::chrono::duration<double>(5.0)) {
	std::cout << "Collide::collide_thread_fork: T=" << tnow << " ratio="
		  << wctduration.count()/cpu_duration << std::endl;
      }
    } else {
      std::cout << "Collide::collide_thread_fork: T=" << tnow << " "
		<< cpu_duration << " CPU sec, "
		<< wctduration.count() << " Wall sec, "
		<< wctduration.count()/cpu_duration
		<< " ratio" << std::endl;
    }
  }

}


double   Collide::EPSMratio = -1.0;
unsigned Collide::EPSMmin   = 0;

std::vector<double> Collide::atomic_weights;

//! Weights in atomic mass units
void Collide::atomic_weights_init()
{
  atomic_weights.resize(15, -1.0);

  atomic_weights[0]  = 0.000548579909; // Mass of electron
  atomic_weights[1]  = 1.0079;	       // Hydrogen
  atomic_weights[2]  = 4.0026;	       // Helium
  atomic_weights[3]  = 6.941;	       // Lithum
  atomic_weights[4]  = 9.0122;	       // Beryllium
  atomic_weights[5]  = 10.811;	       // Boron
  atomic_weights[6]  = 12.011;	       // Carbon
  atomic_weights[7]  = 14.007;	       // Nitrogen
  atomic_weights[8]  = 15.999;	       // Oxygen
  atomic_weights[9]  = 18.998;	       // Florine
  atomic_weights[10] = 20.180;	       // Neon
  atomic_weights[11] = 22.990;	       // Sodium
  atomic_weights[12] = 24.305;	       // Magnesium
  atomic_weights[13] = 26.982;	       // Aluminium
  atomic_weights[14] = 28.085;	       // Silicon
}  

Collide::Collide(ExternalForce *force, Component *comp,
		 double hDiam, double sCross,
		 const std::string& name_id,
		 const std::string& version_id,
		 int nth)
{
  caller = force;
  c0     = comp;
  nthrds = nth;
  
  // For threading
  //
  tlock = PTHREAD_MUTEX_INITIALIZER;

  // Cache the calling tree
  tree = c0->Tree();

  // Initialize atomic weights map
  if (atomic_weights.size()==0) atomic_weights_init();

  // Unintialized molecular weight
  mol_weight = -1.0;

  // Counts the total number of collisions
  colcntT.resize(nthrds);
  
  // Total number of particles processsed
  numcntT.resize(nthrds);
  
  // Total velocity dispersion (i.e. mean temperature)
  tdispT.resize(nthrds);
  
  // Make return cross-section map
  retCrs.resize(nthrds);

  // Number of collisions with inconsistencies (only meaningful for LTE)
  error1T.resize(nthrds, 0);
  
  // Number of particles selected for collision
  sel1T.resize(nthrds, 0);
  
  // Number of particles actually collided
  col1T.resize(nthrds, 0);
  
  // Number of particles processed by the EPSM algorithm
  epsm1T.resize(nthrds, 0);
  
  // Number of cells processed by the EPSM algorithm
  Nepsm1T.resize(nthrds, 0);
  
  // Total mass of processed particles
  tmassT.resize(nthrds, 0);
  
  // True energy lost to dissipation (i.e. radiation)
  decolT.resize(nthrds, 0);
  
  // Full energy lost to dissipation (i.e. radiation) 
  decelT.resize(nthrds, 0);
  
  // Energy excess (true energy)
  exesCT.resize(nthrds, 0);
  
  // Energy excess (full energy)
  exesET.resize(nthrds, 0);
  
  // NTC statistics
  ntcOvr.resize(nthrds, 0);
  ntcAcc.resize(nthrds, 0);
  ntcTot.resize(nthrds, 0);
  ntcVal.resize(nthrds);
  wgtVal.resize(nthrds);
  for (auto &v : ntcVal) v = std::make_shared<circBuf>(bufCap);
  for (auto &v : wgtVal) v = std::make_shared<circBuf>(bufCap);

  if (MFPDIAG) {
    // List of ratios of free-flight length to cell size
    tsratT.resize(nthrds);
    
    // List of fractional changes in KE per cell
    keratT.resize(nthrds);
    
    // List of cooling excess to KE per cell
    deratT.resize(nthrds);
    
    // List of densities in each cell
    tdensT.resize(nthrds);
    
    // List of cell volumes
    tvolcT.resize(nthrds);
    
    // Temperature per cell; assigned in derived class instance
    ttempT.resize(nthrds);
    
    // List of change in energy per cell due to cooling (for LTE only)
    tdeltT.resize(nthrds);
    
    // List of collision selections per particle
    tselnT.resize(nthrds);
    
    // List of cell diagnostic info per cell
    tphaseT.resize(nthrds);
    
    // List of mean-free path info per cell
    tmfpstT.resize(nthrds);
  }
  
  cellist.resize(nthrds);
  
  hsdiam    = hDiam;
  crossfac  = sCross;
  
  seltot    = 0;	      // Count estimated collision targets
  coltot    = 0;	      // Count total collisions
  errtot    = 0;	      // Count errors in inelastic computation
  epsmcells = 0;	      // Count cells in EPSM regime
  epsmtot   = 0;	      // Count particles in EPSM regime
  
  // EPSM diagnostics
  lostSoFar_EPSM = vector<double>(nthrds, 0.0);
  
  if (EPSMratio> 0) use_epsm = true;
  else              use_epsm = false;
  
  stepcount = 0;
  bodycount = 0;
  
  listTime   .resize(nthrds);
  initTime   .resize(nthrds);
  collTime   .resize(nthrds);
  elasTime   .resize(nthrds);
  stat1Time  .resize(nthrds);
  stat2Time  .resize(nthrds);
  stat3Time  .resize(nthrds);
  coolTime   .resize(nthrds);
  cellTime   .resize(nthrds);
  curcTime   .resize(nthrds);
  epsmTime   .resize(nthrds);
  listSoFar  .resize(nthrds);
  initSoFar  .resize(nthrds);
  collSoFar  .resize(nthrds);
  elasSoFar  .resize(nthrds);
  cellSoFar  .resize(nthrds);
  curcSoFar  .resize(nthrds);
  epsmSoFar  .resize(nthrds);
  stat1SoFar .resize(nthrds);
  stat2SoFar .resize(nthrds);
  stat3SoFar .resize(nthrds);
  coolSoFar  .resize(nthrds);
  collCnt    .resize(nthrds, 0);
  
  EPSMT      .resize(nthrds);
  EPSMTSoFar .resize(nthrds);
  for (int n=0; n<nthrds; n++) {
    EPSMT[n].resize(nEPSMT);
    EPSMTSoFar[n].resize(nEPSMT);
  }  
  
  if (TSDIAG) {
    // Accumulate distribution log ratio of flight time to time step
    tdiag  .resize(numdiag, 0);
    tdiag1 .resize(numdiag, 0);
    tdiag0 .resize(numdiag, 0);
    tdiagT .resize(nthrds);
    
    // Accumulate distribution log energy overruns
    Eover  .resize(numdiag, 0);
    Eover1 .resize(numdiag, 0);
    Eover0 .resize(numdiag, 0);
    EoverT .resize(nthrds);
  }
  
  // Accumulate the ratio cooling time to time step each cell
  tcool  .resize(numdiag, 0);
  tcool1 .resize(numdiag, 0);
  tcool0 .resize(numdiag, 0);
  tcoolT .resize(nthrds);
  
  if (MFPCL) mfpCLdata = std::vector< std::vector<double> > (nthrds);

  if (VOLDIAG) {
    Vcnt  .resize(nbits, 0);
    Vcnt1 .resize(nbits, 0);
    Vcnt0 .resize(nbits, 0);
    VcntT .resize(nthrds);
    Vdbl  .resize(nbits*nvold, 0.0);
    Vdbl1 .resize(nbits*nvold, 0.0);
    Vdbl0 .resize(nbits*nvold, 0.0);
    VdblT .resize(nthrds);
  }
  
  for (int n=0; n<nthrds; n++) {
    if (TSDIAG) {
      tdiagT[n] .resize(numdiag, 0);
      EoverT[n] .resize(numdiag, 0);
    }
    if (VOLDIAG) {
      VcntT[n] .resize(nbits, 0);
      VdblT[n] .resize(nbits*nvold, 0.0);
    }
    tcoolT[n] .resize(numdiag, 0);
    tdispT[n] .resize(3, 0);
  }
  
  disptot .resize(3, 0);
  masstot = 0.0;

  use_Eint = -1;
  use_temp = -1;
  use_dens = -1;
  use_delt = -1;
  use_exes = -1;
  use_Kn   = -1;
  use_St   = -1;
  
  if (MFPDIAG) {
    prec .resize(nthrds);
    for (int n=0; n<nthrds; n++)
      prec[n].second .resize(Nmfp, 0);
  }
  
  if (VERBOSE>5) {
    timer_list .resize(2*nthrds);
  }
  
  forkSum  .resize(3);
  snglSum  .resize(3);
  waitSum  .resize(3);
  diagSum  .resize(3);
  joinSum  .resize(3);
  
  if (TIMING) {
    listSum  .resize(3);
    initSum  .resize(3);
    collSum  .resize(3);
    elasSum  .resize(3);
    cellSum  .resize(3);
    epsmSum  .resize(3);
    stat1Sum .resize(3);
    stat2Sum .resize(3);
    stat3Sum .resize(3);
    coolSum  .resize(3);
    numbSum  .resize(3);
  }
  
  EPSMtime .resize(nEPSMT);
  CPUH .resize(12);
  
  // Debug maximum work per cell
  //
  minUsage .resize(nthrds*2, std::numeric_limits<long>::max());
  maxUsage .resize(nthrds*2, 0);
  minPart  .resize(nthrds*2, -1);
  maxPart  .resize(nthrds*2, -1);
  minCollP .resize(nthrds*2, -1);
  maxCollP .resize(nthrds*2, -1);
  
  effortAccum  = false;
  effortNumber .resize(nthrds);

  // Set the threshold from the default value
  //
  use_ntcdb = true;
  ntcThresh = ntcThreshDef;
  ntcFactor = 1.0;

  // Log file identification info
  //
  if (myid==0) {
    std::cout << printDivider << std::endl
	      << "** Collide routine "
	      << name_id << ", version " << version_id << std::endl
	      << printDivider << std::endl << std::endl
	      << printDivider << std::endl
	      << "--- Tree volume = " << tree->Volume() << std::endl
	      << printDivider << std::endl << std::endl;
  }

  // Initialize diagnostic counters
  // 
  pre_collide_diag();
}

Collide::~Collide()
{
  // NADA
}

void Collide::debug_list()
{
  unsigned ncells = tree->Number();
  pHOT_iterator c(*tree);

  for (int cid=0; cid<numprocs; cid++) {

    // Walk through processes
    //
    if (myid == cid) {

      unsigned cellC = 0;
      for (auto l : cellist) cellC += l.size();

      std::ostringstream sout;
      sout << "==== Collide [" << myid << "] level=" << mlev
	   << " active cells/total=" << cellC << "/" << ncells << " ";

      std::cout << std::setw(70)     << std::setfill('=') << std::left
		<< sout.str()        << std::endl
		<< std::setfill(' ') << std::endl;

      if (cellC > 0) {

	//
	//  +--- Summary info
	//  |
	//  V
	if (false) {
	  std::cout << std::setw( 8) << std::right << "TID"
		    << std::setw(12) << std::right << "# cells"
		    << std::setw(12) << std::right << "# bodies"
		    << std::endl
		    << std::setw( 8) << std::right << "---"
		    << std::setw(12) << std::right << "--------"
		    << std::setw(12) << std::right << "--------"
		    << std::endl;
	  
	  for (int n=0; n<nthrds; n++) {
	    unsigned bodCount = 0;
	    for (auto c : cellist[n]) bodCount += c->bods.size();
	    cout << setw(8)  << n
		 << setw(12) << cellist[n].size()
		 << setw(12) << bodCount
		 << std::endl;
	  }
	  
	}      // END: summary info

	else { // BEG: per cell info

	  
	  for (int n=0; n<nthrds; n++) {
	    int nbeg = ncells*(n  )/nthrds;
	    int nend = ncells*(n+1)/nthrds;
	    for (int j=nbeg; j<nend; j++) {
	      int tnum = c.nextCell();
	      std::cout << std::setw(8)  << j
			<< std::setw(12) << c.Cell()
			<< std::setw(12) << cellist[n][j-nbeg]
			<< std::setw(12) << c.Cell()->bods.size()
			<< std::setw(12) << tnum << endl;
	    }

	  } // END: thread loop

	} // END: per cell info

      } // END: cellC>0

    } // process selection

    MPI_Barrier(MPI_COMM_WORLD);

  } // END::process loop
}


const Collide::UU&
Collide::collide(pHOT& tree, sKeyDmap& Fn, int mlevel, bool diag)
{
  std::get<0>(ret) = std::get<1>(ret) = 0;
  mlev = mlevel;

  snglTime.start();

  // Make cellist
  // 
  for (int n=0; n<nthrds; n++) cellist[n].clear();
  ncells = 0;
  
  // For debugging
  //
  unsigned nullcell = 0, totalcell = 0, totalbods = 0;
  
  for (unsigned M=mlevel; M<=multistep; M++) {
    
    // Don't queue null cells
    //
    if (tree.CLevels(M).size()) {

      for (auto ic : tree.CLevels(M)) {
	if (ic->bods.size()) {
	  cellist[(ncells++)%nthrds].push_back(ic);
	  bodycount += ic->bods.size();
	  totalbods += ic->bods.size();
	} else {
	  nullcell++;
	}
	totalcell++;
      }
    }
  }
  stepcount++;
  
  unsigned allbods;
  MPI_Allreduce(&totalbods, &allbods, 1, MPI_UNSIGNED, MPI_SUM, 
		MPI_COMM_WORLD);

  // ret<0> Will contain sum of all particles at this level
  //
  std::get<0>(ret) = allbods;

  {
    std::ostringstream sout;
    sout << "Collide::collide: AFTER totalbods allreduce, N="
	 << std::get<0>(ret) << ", Ntest=" << allbods;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }


  if (CDEBUG) {

    if (myid==0) {
      std::cout << std::endl << "***" << std::endl << std::left 
		<< "*** T = "  << std::setw(10) << tnow 
		<< " Level = " << std::setw(10) << mlevel
		<< " Nbods = " << std::setw(10) << std::get<0>(ret)
		<< " dT = "    << dtime << std::endl
		<< "***" << std::endl;

      std::cout << std::endl << printDivider << std::endl
		<< std::left 
		<< std::setw(8)  << "Node"
		<< std::setw(6)  << "Level"
		<< std::setw(10) << "Bodies"
		<< std::setw(10) << "dT"
		<< std::endl << printDivider << std::endl;
    }
		 
    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	std::cout << "#" << std::setw(4) << std::left << myid << " | "
		  << std::setw(6)  << mlevel 
		  << std::setw(10) << totalbods
		  << std::setw(10) << dtime
		  << std::endl << std::right;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if (myid==0) std::cout << printDivider << std::endl;
  }

  if (CDEBUG) {
    if (nullcell)
      std::cout << "DEBUG: null cells " << nullcell << "/" 
		<< totalcell << std::endl;
    
    debug_list();
  }
  snglTime.stop();
  
  // Needed for meaningful timing results
  //
  waitTime.start();
  (*barrier)("Collide::collide: in sync BEFORE collision fork",
	     __FILE__, __LINE__);
  MPI_Barrier(MPI_COMM_WORLD);
  waitTime.stop();
  
  // For effort debugging
  //
  if (mlevel==0) effortAccum = true;
  
  forkTime.start();
  if (0 && totalbods) {
    ostringstream sout;
    sout << "before fork, " << __FILE__ << ": " << __LINE__;
    tree.checkBounds(2.0, sout.str().c_str());
  }

  collide_thread_fork(&Fn);

  if (0 && totalbods) {
    ostringstream sout;
    sout << "after fork, " << __FILE__ << ": " << __LINE__;
    tree.checkBounds(2.0, sout.str().c_str());
  }
  forkSoFar = forkTime.stop();
  
  snglTime.start();
  if (diag) {
				// Collect
    std::get<1>(ret) = post_collide_diag();
    pre_collide_diag();		// Reinitialize
  }
  snglSoFar = snglTime.stop();
  
  (*barrier)("Collide::collide: AFTER collision fork",
	     __FILE__, __LINE__);

  // Effort diagnostics
  //
  if (mlevel==0 && EFFORT && effortAccum) {

    (*barrier)("Collide::collide: BEFORE effort diagnostic",  
	       __FILE__, __LINE__);

    static double lastTime = -1.0;
    static double timeIntv =  1.0e-10;

    unsigned char doweet = 0;
    if (tnow - lastTime > timeIntv) doweet = 1;
    MPI_Bcast(&doweet, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

    if (doweet) {

      (*barrier)("Collide::collide: DOING effort diagnostic",  
		 __FILE__, __LINE__);

      // Compute summary statistics
      //
      double mean=0.0, var2=0.0, cnts=0.0;
      for (int n=0; n<nthrds; n++) {
	for (auto it : effortNumber[n]) {
	  double val = static_cast<double>(it.first)/it.second;
	  mean += val;
	  var2 += val*val;
	  cnts += 1;
	}
      }

      if (cnts>0.0) {
	mean /= cnts;
	var2  = var2/cnts - mean*mean;
      }
      
      // Send into to root
      //
      std::vector<double> mean_all(numprocs);
      MPI_Gather(&mean, 1, MPI_DOUBLE, &mean_all[0], 1, MPI_DOUBLE, 0,
		 MPI_COMM_WORLD);

      std::vector<double> var2_all(numprocs);
      MPI_Gather(&var2, 1, MPI_DOUBLE, &var2_all[0], 1, MPI_DOUBLE, 0,
		 MPI_COMM_WORLD);
      
      std::vector<double> cnts_all(numprocs);
      MPI_Gather(&cnts, 1, MPI_DOUBLE, &cnts_all[0], 1, MPI_DOUBLE, 0,
		 MPI_COMM_WORLD);
      
      (*barrier)("Collide::collide: effort diagnostic GATHER complete",  
		 __FILE__, __LINE__);

      //
      // Root nodes writes the file in node order
      //
      if (myid==0) {
	
	std::ostringstream ostr;
	ostr << outdir << runtag << ".collide.effort";

	ofstream out;
	try {
	  out.open(ostr.str().c_str(), ios::app);
	}
	catch (std::ios_base::failure& e) {
	  std::cerr << "Collide::collide, error opening <"
		    << ostr.str() << ">: " << e.what() << std::endl;
	}

	if (out.good()) {
	  static bool firstTime = true;
	  
	  if (firstTime) {
	    out << std::setw( 6) << "Pid"
		<< std::setw(18) << "Time"
		<< std::setw(18) << "Mean(musec)"
		<< std::setw(18) << "Var(musec)"
		<< std::setw(18) << "Counts"
		<< std::endl
		<< std::setw( 6) << "-----"
		<< std::setw(18) << "----------"
		<< std::setw(18) << "----------"
		<< std::setw(18) << "----------"
		<< std::setw(18) << "----------"
		<< std::endl;
	    
	    firstTime = false;
	  }
	  
	  for (int i=0; i<numprocs; i++) {
	    out << std::setw( 6) << i
		<< std::setw(18) << tnow
		<< std::setw(18) << mean_all[i]
		<< std::setw(18) << sqrt(fabs(var2_all[i]))
		<< std::setw(18) << cnts_all[i]
		<< std::endl;
	  }
	  
	} else {
	  cerr << "Process " << myid 
	       << ": error opening <" << ostr.str() << ">" << endl;
	}

	out.close();
      }

      // Reset the list
      //
      effortAccum = false;
      for (int n=0; n<nthrds; n++) effortNumber[n].clear();
    }

    (*barrier)("Collide::collide: AFTER effort diagnostic",  
	       __FILE__, __LINE__);
  }
  
  caller->print_timings("Collide: collision thread timings", timer_list);
  
  // Persist NTC database
  //
  if (mlevel==0 and use_ntcdb) ntcdb.update();

  {
    std::ostringstream sout;
    sout << "Collide::collide: EXITING, bods=" << std::get<0>(ret)
	 << ", allbods=" << allbods;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

  return ret;
}

void Collide::dispersion(vector<double>& disp)
{
  disp = disptot;
  if (masstot>0.0) {
    for (unsigned k=0; k<3; k++) disp[k] /= masstot;
  }
  for (unsigned k=0; k<3; k++) disptot[k] = 0.0;
  masstot = 0.0;
}


void * Collide::collide_thread(void * arg)
{
  sKeyDmap *Fn  = static_cast<sKeyDmap*>(((thrd_pass_arguments*)arg)->fn  );
  int id        = static_cast<int>      (((thrd_pass_arguments*)arg)->id  );
  
  thread_timing_beg(id);
  
  if (id==0) {
    std::ostringstream sout;
    sout << "Collide::collide: ENTERING collide_thread, T=" << tnow;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

  // Initialize cell loop diagnostics
  //
  pre_cell_loop(id);

  // Start execution timer
  //
  cellTime[id].start();
  
  // DEEP DEBUG
  if (false) {
    unsigned elem = 0, celltot = 0, cellsum = 0;
    for (auto v : cellist) {
      elem++;
      celltot += v.size();
    }

    std::cout << "[" << myid << "] cells=" << celltot
	      << "/" << elem << std::endl;

    MPI_Reduce(&celltot, &cellsum, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) std::cout << "[sum] cells=" << cellsum << std::endl;
  }
  // END DEBUG

  // Store counters for collisional excitation counting (for TESTING)
  //
#ifdef TESTSPREAD
  static std::tuple<double, double, double, double, int, int, int> running_tally = {0, 0, 0, 0, 0, 0, 0};
  //                                                                                ^  ^  ^  ^  ^  ^  ^
  // 0 Accepted collisions----------------------------------------------------------+  |  |  |  |  |  |
  // 1 Proposed collisions-------------------------------------------------------------+  |  |  |  |  |
  // 2 Cross section----------------------------------------------------------------------+  |  |  |  |
  // 3 Target ratio--------------------------------------------------------------------------+  |  |  |
  // 4 Cumulated pairs--------------------------------------------------------------------------+  |  |
  // 5 Number of trials----------------------------------------------------------------------------+  |
  // 6 Number accepted--------------------------------------------------------------------------------+
#endif

  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {
    
    // The current cell
    //
    pCell *c = cellist[id][j];

    // Skip cell if this time has already been computed
    if (c->time >= tnow) {
      continue;
    }

#ifdef XC_DEEPT
  std::cout << "**TIME=" << tnow << std::endl;
#endif

#ifdef USE_GPTL
    GPTLstart("Collide::bodylist");
#endif
    
    int EPSMused = 0;
    
    // Start the effort time
    //
    curcTime[id].reset();
    curcTime[id].start();
    listTime[id].start();
    
    // Number of particles in this cell
    //
    unsigned number = c->bods.size();
    numcntT[id].push_back(number);
    
#ifdef USE_GPTL
    GPTLstop("Collide::bodylist");
#endif
    
    // Skip cells with only one particle
    //
    if ( number < 2 ) {
      colcntT[id].push_back(0);
      // Stop timers
      curcTime[id].stop();
      listSoFar[id] = listTime[id].stop();
      // Skip to the next cell
      continue;
    }
    
#ifdef USE_GPTL
    GPTLstart("Collide::prelim");
    GPTLstart("Collide::energy");
#endif
    
    listSoFar[id] = listTime[id].stop();
    stat1Time[id].start();
    
    // Energy lost in this cell
    //
    decolT[id] = 0.0;
    decelT[id] = 0.0;
    
    // Compute 1.5 times the mean relative velocity in each MACRO cell
    //
    pCell *samp = c->sample;

    sKeyUmap::iterator it1, it2;

    //
    // Sanity check
    //
    if (samp == 0x0) {
      cout << "Process "  << myid << " in collide: no sample cell"
	   << ", owner="   << c->owner << hex
	   << ", mykey="   << c->mykey
	   << ", mask="    << c->mask  << dec
	   << ", level="   << c->level    
	   << ", Count="   << c->ctotal
	   << ", maxplev=" << c->maxplev;
      if (tree->onFrontier(c->mykey)) cout << ", ON frontier" << endl;
      else cout << ", NOT on frontier" << endl;
      
    }

    double crm = 0.0;
    
    if (SampleCRM and samp) {
      if (samp->stotal[0]>0.0) {
	for (unsigned k=0; k<3; k++) {
	  crm += (samp->stotal[1+k] - 
		  samp->stotal[4+k]*samp->stotal[4+k]/samp->stotal[0])
	    /samp->stotal[0];
	}
      }
      crm  = crm>0.0 ? sqrt(2.0*crm) : 0.0;
    } else {
      if (c->stotal[0]>0.0) {
	for (unsigned k=0; k<3; k++) {
	  crm += (c->stotal[1+k] - 
		  c->stotal[4+k]*c->stotal[4+k]/c->stotal[0])
	    /c->stotal[0];
	}
      }
      crm  = crm>0.0 ? sqrt(2.0*crm) : 0.0;
    }
    
    stat1SoFar[id] = stat1Time[id].stop();
    stat2Time [id].start();
    
#ifdef USE_GPTL
    GPTLstop ("Collide::energy");
    GPTLstart("Collide::mfp");
#endif
    
    // KE in the cell
    //
    double kedsp=0.0;
    if (MFPDIAG) {
      if (c->stotal[0]>0.0) {
	for (unsigned k=0; k<3; k++) 
	  kedsp += 
	    0.5*(c->stotal[1+k] - c->stotal[4+k]*c->stotal[4+k]/c->stotal[0]);
      }
    }
    
    // Timestep for this cell
    //
    double tau = dtime / (1<<c->maxplev);

    // Volume in the cell
    //
    double volc = c->Volume();
    
    // Mass in the cell
    //
    double mass = c->Mass();
    
    // Mass density in the cell
    double dens = mass/volc;
    
    if (mass <= 0.0) continue;
    
    // Cell length
    //
    double cL = pow(volc, 0.33333333);
    
    
    // Cell initialization
    //
    initialize_cell(c, crm, tau, id);

    // Number of bodies in the cell
    //
    int nbods = c->bods.size();

    // Per species quantities
    //
    double         meanLambda, meanCollP, totalNsel;
    NTC::InteractD crossIJ = totalCrossSections(c, crm, id);
    NTC::InteractD nselM   = generateSelection (c, crm, Fn, tau,
						meanLambda, meanCollP, totalNsel, id);
    

    // True for intensive debugging only
    //  +------------------------------+
    //  |
    //  v
    if (false) {
      
      const double cunit = 1e-14/(TreeDSMC::Lunit*TreeDSMC::Lunit);

      // nselM is an NTC::InteractD. NTC::InteractD is an interface to
      // a std::map<NTC::T, double>.  NTC::T is a tuple of the
      // interaction type and the two speciesKey values.

      std::cout << std::left << setw(20) << "Species"
		<< std::left << setw(30) << "Interaction"
		<< std::left << setw(18) << "Cross section"
		<< std::left << setw(18) << "Selection"
		<< std::endl
		<< std::left << setw(20) << "-------"
		<< std::left << setw(30) << "-----------"
		<< std::left << setw(18) << "-------------"
		<< std::left << setw(18) << "---------"
		<< std::endl;

      double totalNsel = 0.0;
      for (auto v : nselM.v) {
	std::ostringstream sout1;
	sout1 << "[(" << (short)std::get<1>(v.first).first
	      << ", " << (short)std::get<1>(v.first).second
	      << ")(" << (short)std::get<2>(v.first).first
	      << ", " << (short)std::get<2>(v.first).second
	      << ")]";

	std::cout << std::left << std::setw(20) << sout1.str()
		  << std::left << std::setw(30) << labels[std::get<0>(v.first)]
		  << std::left << std::setw(18) << crossIJ[v.first]()
		  << std::left << std::setw(18) << v.second()
		  << std::endl;

	totalNsel += v.second();
      }
      std::cout << "Total Nsel=" << totalNsel << std::endl;
    }

#ifdef USE_GPTL
    GPTLstop ("Collide::mfp");
#endif
    
#ifdef USE_GPTL
    GPTLstart("Collide::mfp_diag");
#endif
    
    if (MFPDIAG) {
      
      // Diagnostics
      //
      tsratT[id].push_back(crm*tau/pow(volc,0.33333333));
      tdensT[id].push_back(dens);
      tvolcT[id].push_back(volc);
      
      double posx, posy, posz;
      c->MeanPos(posx, posy, posz);
      
      // MFP/side = MFP/vol^(1/3)
      
      prec[id].first = meanLambda/pow(volc, 0.33333333);
      prec[id].second[0] = sqrt(posx*posx + posy*posy);
      prec[id].second[1] = posz;
      prec[id].second[2] = sqrt(posx*posx+posy*posy*+posz*posz);
      prec[id].second[3] = mass/volc;
      prec[id].second[4] = volc;
      
      tmfpstT[id].push_back(prec[id]);
    }
    
    // Ratio of cell size to mean free path
    //
    double mfpCL = cL/meanLambda;

    if (MFPCL) mfpCLdata[id].push_back(mfpCL);
    
    if (TSDIAG) {		// Diagnose time step in this cell
      double vmass;
      std::vector<double> V1, V2;
      c->Vel(vmass, V1, V2);
      double scale = c->Scale();
      double taudiag = 1.0e40;
      for (int k=0; k<3; k++) {	// Time of flight
	taudiag = min<double>
	  (pHOT::sides[k]*scale/(fabs(V1[k]/vmass)+sqrt(V2[k]/vmass)+1.0e-40), 
	   taudiag);
      }
      
      int indx = (int)floor(log(taudiag/tau)/log(4.0) + 5);
      if (indx<0 ) indx = 0;
      if (indx>10) indx = 10;
      tdiagT[id][indx]++;
    }
    
    if (VOLDIAG) {
      if (c->level<nbits) {
	VcntT[id][c->level]++;
	VdblT[id][c->level*nvold+0] += dens;
	VdblT[id][c->level*nvold+1] += 1.0 / mfpCL;
	VdblT[id][c->level*nvold+2] += meanCollP;
	VdblT[id][c->level*nvold+3] += crm*tau / cL;
	VdblT[id][c->level*nvold+4] += number;
	VdblT[id][c->level*nvold+5] += number*number;
      }
    }

    stat2SoFar[id] = stat2Time[id].stop();
    
    // Species map for collisions

    unsigned colc = 0;
    std::map<speciesKey, std::vector<unsigned long> > bmap;

#ifdef USE_GPTL
    GPTLstop ("Collide::mfp_diag");
    GPTLstart("Collide::cell_init");
#endif
    
    
#ifdef USE_GPTL
    GPTLstop("Collide::cell_init");
    GPTLstop("Collide::prelim");
#endif
    // No collisions, primarily for testing . . .
    //
    if (DRYRUN) continue;
    
    collTime[id].start();
    collCnt[id]++;
    // Number of collisions per particle:
    // assume equipartition if large
    //
    if (use_epsm && meanCollP > EPSMratio && number > EPSMmin) {
      
      EPSMused = 1;
      
#ifdef USE_GPTL
      GPTLstart("Collide::cell_init");
#endif
      initTime[id].start();
      initialize_cell_epsm(c, nselM, crm, tau, id);
      initSoFar[id] = initTime[id].stop();
      
#ifdef USE_GPTL
      GPTLstop("Collide::cell_init");
      GPTLstart("Collide::EPSM");
#endif
      epsmTime[id].start();
      EPSM(c, id);
      epsmSoFar[id] = epsmTime[id].stop();
#ifdef USE_GPTL

      GPTLstop ("Collide::EPSM");
#endif
      
    } else {
      
#ifdef USE_GPTL
      GPTLstart("Collide::cell_init");
#endif
      initTime[id].start();
      initialize_cell_dsmc(c, nselM, crm, tau, id);
      initSoFar[id] = initTime[id].stop();
      
#ifdef USE_GPTL
      GPTLstop("Collide::cell_init");
      GPTLstart("Collide::inelastic");
#endif
      for (size_t k=0; k<c->bods.size(); k++) {
	unsigned long kk = c->bods[k];
	Particle* p = tree->Body(kk);

	speciesKey skey = p->skey;
	if (use_key>=0 and skey==Particle::defaultKey)
	  skey = p->skey = KeyConvert(p->iattrib[use_key]).getKey();

	bmap[skey].push_back(kk);
      }
    }
    
    int acceptCount = 0;
    
    for (auto v : nselM.v) {

      int totalCount = ceil(v.second());
      
      NTC::T T = v.first;

      int        TT = std::get<0>(T); // For debugging only
      speciesKey i1 = std::get<1>(T);
      speciesKey i2 = std::get<2>(T);
	
      sKeyPair k(i1, i2);

#ifdef XC_DEEP12
      printf("NPAIR=%8d NSEL=%13.6e T=%d\n", totalCount, v.second(), TT);
#endif

#ifdef TESTSPREAD
      double tfrac = 1.0;
      if (TT>=6 and v.second()<TestSpreadCount) {
	totalCount = TestSpreadCount;
	tfrac = v.second()/totalCount;
      }

      if (TT == 7) {
	std::get<1>(running_tally) += v.second();
	std::get<4>(running_tally) += nbods*(nbods-1)/2;
      }
#endif

      for (int np=0; np<totalCount; np++) {

	// Fractional check
	//
	double frc = v.second() - np, wght = 1.0;
	double  R0 = unit(random_gen);
#ifdef  XC_DEEP12
	printf("FRC=%13.6e R=%13.6e T=%d\n", frc, R0, TT);
#endif

#ifdef  TESTSPREAD
	wght = tfrac;
#else

#ifdef  WEIGHTED
	if (frc < 1.0) wght = frc;
#else
	if (frc < 1.0 and R0 > frc) break;
	//      ^            ^
	//      |            |
	//      |            +--- select with probability frc
	//      |
	//      +--- Only use fractional part on final candidate
#endif

#endif

	// Pick a pair of particles from the cell
	//
	auto pr = pairSelect(id, T, nbods);

	int i1 = pr.first;
	int i2 = pr.second;

	// Get index from body map for the cell
	//
	Particle* const p1 = tree->Body(c->bods[i1]);
	Particle* const p2 = tree->Body(c->bods[i2]);
	
	// Calculate pair's relative speed (pre-collision)
	//
	vector<double> crel(3);
	double cr = 0.0;
	for (int j=0; j<3; j++) {
	  crel[j] = p1->vel[j] - p2->vel[j];
	  cr += crel[j]*crel[j];
	}
	cr = sqrt(cr);

	// No point in inelastic collsion for zero velocity . . . 
	//
	if (cr == 0.0) continue;

	// Per particle initialization for cross section
	//
	pairInfo(id, c, p1, p2, cr);
	
	// Accept or reject candidate pair according to relative speed
	//
	const double cunit = 1e-14/(TreeDSMC::Lunit*TreeDSMC::Lunit);
	double Cross = crossSection(id, c, p1, p2, cr, T);
	bool ok = false;
	  
	double crsvel  = crossIJ[v.first]()/cunit;
	if (use_ntcdb and ntcdb[samp->mykey].Ready(k, T)) {
	  crsvel = ntcdb[samp->mykey].CrsVel(k, T, ntcThresh) * ntcFactor;
	}
	
	if (crsvel > 0.0) {

	  double scrs = Cross / cunit;
	  double prod = cr * scrs;
	  double targ = prod / crsvel;
	  
#ifdef XC_DEEP12
	  std::cout << "TARG=" << targ << " T=" << TT << std::endl;
#endif
	  if (NTC or NTCnodb)
	    ok = ( targ > unit(random_gen) );
	  else
	    ok = true;

#ifdef TESTSPREAD
	  if (TT == 7 and Cross > 0.0) {
	    std::get<0>(running_tally) += tfrac;
	    std::get<2>(running_tally) += Cross;
	    std::get<3>(running_tally) += targ;
	    std::get<5>(running_tally) += 1;
	    if (ok) std::get<6>(running_tally) += 1;

	  }
#endif


#ifdef XC_DEEP14
	  static long tally = 0, accept = 0;
	  if (TT==7) {
	    if (ok) accept++; tally++;
	    std::cout << "col_excite: targ=" << targ
		      << " accept=" << accept
		      << " tally=" << tally
		      << " ratio=" << (double)accept/tally
		      << std::endl;
	  }
#endif
	  // Update v_max and cross_max for NTC
	  //
	  if (NTC) {
				// Diagnostic
	    ntcVal[id]->push_back(targ);
				// Over NTC max average
	    if (targ >= 1.0) ntcOvr[id]++;
				// Used / Total
	    if (ok) ntcAcc[id]++;
	    ntcTot[id]++;

				// Accumulate average
	    pthread_mutex_lock(&tlock);
	    ntcdb[samp->mykey].Add(k, T, prod);
	    pthread_mutex_unlock(&tlock);
	    if (T == NTC::single) {
	      std::cout << "ntc singleton" << std::endl;
	    }

				// Sanity check
	    if (numSanityMsg and 
		ntcAcc[id] >= numSanityVal and 	
		ntcAcc[id] %  numSanityFreq == 0) {
	      
	      std::ostringstream sout;
	      sout << "<"
		   << k.first .first << "," << k.first .second << "|"
		   << k.second.first << "," << k.second.second << ">";
	      
	      std::cout << "Proc " << myid << " thread=" << id 
			<< ": cell="  << c->mykey
			<< ", count=" << c->bods.size()
			<< ", targ="  << targ
			<< ", mfpCL=" << mfpCL
			<< ", nselM=" << v.second()
			<< ", ntcOvr=" << ntcOvr[id]
			<< " for "    << sout.str() << std::endl
			<< " has logged " << ntcAcc[id] << " collisions!"
			<< " You may wish to cancel this run and"  << std::endl
			<< " adjust the cell size or particle number." 
			<< std::endl;
	    }
	  } // End: NTC
	} else {
	  ok = false;
	}
	    
	if (ok) {

#ifdef XC_DEEP12
	  printf("SELECTED %d\n", TT);
#endif
	  // Counter for time-step selection
	  //
	  acceptCount++;

	  elasTime[id].start();
	    
	  // If pair accepted, select post-collision velocities
	  //
	  colc++;			// Collision counter
	  
	  // Do inelastic stuff
	  //
	  error1T[id] += inelastic(id, c, p1, p2, &cr, wght, T);
	  
	  // Update the particle velocity
	  //
	  velocityUpdate(p1, p2, cr);
	  
	} // Inelastic computation

      } // Particle loop

    } // Species loop
    
    elasSoFar[id] = elasTime[id].stop();
    

    // Count collisions
    //
    colcntT[id].push_back(colc);
    sel1T[id] += acceptCount;
    col1T[id] += colc;
      
    // Selection number per particle
    //
    if (MFPDIAG and nbods>0)
      tselnT[id].push_back(acceptCount/nbods);
    
    // Cache acceptance fraction for scaling MFP for time step
    // selection
    //
    if (acceptCount > 0) {
      /*
      double scale = static_cast<double>(totalCount)/acceptCount; 
      meanLambda *= scale;
      meanCollP  /= scale;
      */
      if (MFPTS) {
	pthread_mutex_lock(&tlock);
	selMFP[c] = meanLambda;
	pthread_mutex_unlock(&tlock);
      }
    }

#ifdef USE_GPTL
    GPTLstop("Collide::inelastic");
#endif
    collSoFar[id] = collTime[id].stop();
  
#ifdef USE_GPTL
    GPTLstart("Collide::diag");
#endif
  
    // Update collision counter in the cell
    //
    c->iattrib["collCount"] = acceptCount;

    // Compute dispersion diagnostics
    //
    stat3Time[id].start();
  
    double tmass = 0.0;
    vector<double> velm(3, 0.0), velm2(3, 0.0);
    for (unsigned j=0; j<number; j++) {
      Particle* p = tree->Body(c->bods[j]);
      for (unsigned k=0; k<3; k++) {
	velm[k]  += p->mass*p->vel[k];
	velm2[k] += p->mass*p->vel[k]*p->vel[k];
      }
      tmass += p->mass;
    }
  
    if (tmass>0.0) {
      for (unsigned k=0; k<3; k++) {
	
	velm[k] /= tmass;
	velm2[k] = velm2[k] - velm[k]*velm[k]*tmass;
	if (velm2[k]>0.0) {
	  tdispT[id][k] += velm2[k];
	  tmassT[id]    += tmass;
	}
      }
    }
    
    //
    // General hook for the derived classes for specific computations
    // and diagnostics
    //
  
    finalize_cell(c, Fn, kedsp, tau, id);
  
    // Update cell time
    //
    c->time = tnow;

    stat3SoFar[id] = stat3Time[id].stop();
  
    //
    // Compute Knudsen and/or Strouhal number
    //
    if (use_Kn>=0 || use_St>=0) {
      double cL = pow(volc, 0.33333333);
      double Kn = meanLambda/cL;
      double St = cL/fabs(tau*sqrt(fabs(kedsp))+1.0e-18);
      for (unsigned j=0; j<number; j++) {
	Particle* p = tree->Body(c->bods[j]);
	if (use_Kn>=0) p->dattrib[use_Kn] = Kn;
	if (use_St>=0) p->dattrib[use_St] = St;
      }
    }
    
#ifdef USE_GPTL
    GPTLstop("Collide::diag");
#endif
  
    // Record effort per particle in microseconds
    //
    curcSoFar[id] = curcTime[id].stop();
    double tt = curcSoFar[id];
    if (EFFORT) {
      if (effortAccum) 
	effortNumber[id].push_back(pair<long, unsigned>(tt, number));
      double effort = static_cast<double>(tt)/number;
      for (unsigned k=0; k<number; k++) 
	tree->Body(c->bods[k])->effort += effort;
    }
  
    // Usage debuging
    //
    if (minUsage[id*2+EPSMused] > tt) {
      minUsage[id*2+EPSMused] = tt;
      minPart [id*2+EPSMused] = number;
      minCollP[id*2+EPSMused] = meanCollP;
    }
    if (maxUsage[id*2+EPSMused] < tt) {
      maxUsage[id*2+EPSMused] = tt;
      maxPart [id*2+EPSMused] = number;
      maxCollP[id*2+EPSMused] = meanCollP;
    }
    
  } // Loop over cells

#ifdef TESTSPREAD
  if (std::get<1>(running_tally)>0.0) {
    std::cout << "COLEXCITE tally:"
	      << " acp:" << std::setw(16) << std::get<0>(running_tally)
	      << " pro:" << std::setw(16) << std::get<1>(running_tally)
	      << " crs:" << std::setw(10) << std::get<4>(running_tally)
	      << " a/p:" << std::setw(16) << std::get<0>(running_tally)/std::get<1>(running_tally)
	      << " cvg: " << std::setw(16) << std::get<2>(running_tally)/std::get<1>(running_tally)
	      << " num: " << std::setw(10) << std::get<6>(running_tally);
    if (std::get<5>(running_tally)>0)
      std::cout << " trg: "
		<< std::setw(16) << std::get<3>(running_tally)/std::get<5>(running_tally);
    std::cout << std::endl;
  }
#endif

  if (id==0) {
    std::ostringstream sout;
    sout << "Collide::collide: AFTER cell loop, T=" << tnow;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

  cellSoFar[id] = cellTime[id].stop();

  // Diagnostics at end of cell loop
  //
  post_cell_loop(id);

  thread_timing_end(id);
  
  return (NULL);
}


unsigned Collide::medianNumber() 
{
  MPI_Status s;
  
  if (myid==0) {
    unsigned num;
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<unsigned> tmp(num);
      MPI_Recv(&tmp[0], num, MPI_UNSIGNED, n, 40, MPI_COMM_WORLD, &s);
      numcnt.insert(numcnt.end(), tmp.begin(), tmp.end());
    }
    
    std::sort(numcnt.begin(), numcnt.end()); 
    
    if (EXTRA) {
      string file = outdir + "tmp.numcnt";
      ofstream out(file.c_str());
      for (unsigned j=0; j<numcnt.size(); j++)
	out << setw(8) << j << setw(18) << numcnt[j] << endl;
    }


    
    if (numcnt.size()) 
      return numcnt[numcnt.size()/2]; 
    else
      return 0;
    
  } else {
    unsigned num = numcnt.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&numcnt[0], num, MPI_UNSIGNED, 0, 40, MPI_COMM_WORLD);
    
    return 0;
  }
}

unsigned Collide::medianColl() 
{ 
  MPI_Status s;
  
  if (myid==0) {
    unsigned num;
    vector<unsigned> coltmp(colcnt);
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<unsigned> tmp(num);
      MPI_Recv(&tmp[0], num, MPI_UNSIGNED, n, 40, MPI_COMM_WORLD, &s);
      coltmp.insert(coltmp.end(), tmp.begin(), tmp.end());
    }
    
    std::sort(coltmp.begin(), coltmp.end()); 
    
    if (EXTRA) {
      ostringstream ostr;
      ostr << outdir << runtag << ".colcnt";
      ofstream out(ostr.str().c_str());
      for (unsigned j=0; j<coltmp.size(); j++)
	out << setw(8) << j << setw(18) << coltmp[j] << endl;
    }
    
    return coltmp[coltmp.size()/2]; 
    
  } else {
    unsigned num = colcnt.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&colcnt[0], num, MPI_UNSIGNED, 0, 40, MPI_COMM_WORLD);
    return 0;
  }
  
}

void Collide::collQuantile(vector<double>& quantiles, vector<double>& coll_)
{
  MPI_Status s;
  
  if (myid==0) {
    unsigned num;
    vector<unsigned> coltmp(colcnt);
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<unsigned> tmp(num);
      MPI_Recv(&tmp[0], num, MPI_UNSIGNED, n, 40, MPI_COMM_WORLD, &s);
      coltmp.insert(coltmp.end(), tmp.begin(), tmp.end());
    }
    
    std::sort(coltmp.begin(), coltmp.end()); 
    
    coll_ .resize(quantiles.size(), 0);
    if (coltmp.size()) {
      for (unsigned j=0; j<quantiles.size(); j++)
	coll_[j] = coltmp[Qi(quantiles[j],coltmp.size())];
    }
    
    ostringstream ostr;
    ostr << outdir << runtag << ".coll_counts";
    ifstream in(ostr.str().c_str());
    in.close();
    if (in.fail()) {
      ofstream out(ostr.str().c_str());
      out << left
	  << setw(14) << "# Time" 
	  << setw(10) << "Quantiles"
	  << setw(10) << "Counts"
	  << endl;
    }
    
    ofstream out(ostr.str().c_str(), ios::app);
    for (unsigned j=0; j<quantiles.size(); j++) {
      out << setw(14) << tnow
	  << setw(10) << quantiles[j] 
	  << setw(10) << coll_[j] << endl;
    }
    out << endl;
    
  } else {
    unsigned num = colcnt.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&colcnt[0], num, MPI_UNSIGNED, 0, 40, MPI_COMM_WORLD);
  }
}

void Collide::mfpsizeQuantile(vector<double>& quantiles, 
			      vector<double>& mfp_, 
			      vector<double>& ts_,
			      vector<double>& coll_,
			      vector<double>& cool_,
			      vector<double>& rate_,
			      unsigned &collnum, unsigned &coolnum) 
{
  if (!MFPDIAG) return;
  
  MPI_Status s;
  
  if (myid==0) {
    unsigned nmt, nmb, num;
    
    // Temporaries so we don't touch the
    // root node's counters
    
    std::vector<Precord> phsI(tphase), mfpI(tmfpst);
    std::vector<double>  ratI(tsrat), denI(tdens), volI(tvolc);
    std::vector<double>  temI(ttemp), selI(tseln), kerI(kerat);
    std::vector<double>  delI(tdelt), derI(derat);
    
    
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&nmt, 1, MPI_UNSIGNED, n, 37, MPI_COMM_WORLD, &s);
      MPI_Recv(&nmb, 1, MPI_UNSIGNED, n, 38, MPI_COMM_WORLD, &s);
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      
      vector<double> tmb(nmb), tmt(nmt), tmp(num);
      
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 40, MPI_COMM_WORLD, &s);
      ratI.insert(ratI.end(), tmp.begin(), tmp.end());
      
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 41, MPI_COMM_WORLD, &s);
      denI.insert(denI.end(), tmp.begin(), tmp.end());
      
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 42, MPI_COMM_WORLD, &s);
      volI.insert(volI.end(), tmp.begin(), tmp.end());
      
      MPI_Recv(&tmp[0], nmt, MPI_DOUBLE, n, 43, MPI_COMM_WORLD, &s);
      temI.insert(temI.end(), tmp.begin(), tmp.end());
      
      MPI_Recv(&tmp[0], nmt, MPI_DOUBLE, n, 44, MPI_COMM_WORLD, &s);
      delI.insert(delI.end(), tmp.begin(), tmp.end());
      
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 45, MPI_COMM_WORLD, &s);
      selI.insert(selI.end(), tmp.begin(), tmp.end());
      
      MPI_Recv(&tmb[0], nmb, MPI_DOUBLE, n, 46, MPI_COMM_WORLD, &s);
      kerI.insert(kerI.end(), tmb.begin(), tmb.end());
      
      MPI_Recv(&tmb[0], nmb, MPI_DOUBLE, n, 47, MPI_COMM_WORLD, &s);
      derI.insert(derI.end(), tmb.begin(), tmb.end());
      
      std::vector<Precord> tmp2(nmt), tmp3(num), phsI(tphase);
      
      MPI_Recv(&tmt[0], nmt, MPI_DOUBLE, n, 48, MPI_COMM_WORLD, &s);
      for (unsigned k=0; k<nmt; k++) {
	// Load density
	tmp2[k].first = tmt[k];
	// Initialize record
	tmp2[k].second .resize(Nphase, 0);
      }
      for (unsigned l=0; l<Nphase; l++) {
	MPI_Recv(&tmp[0], nmt, MPI_DOUBLE, n, 49+l, MPI_COMM_WORLD, &s);
	for (unsigned k=0; k<nmt; k++) tmp2[k].second[l] = tmp[k];
      }
      phsI.insert(phsI.end(), tmp2.begin(), tmp2.end());
      
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 49+Nphase, MPI_COMM_WORLD, &s);
      for (unsigned k=0; k<num; k++) {
	// Load mfp
	tmp3[k].first = tmp[k];
	// Initialize record
	tmp3[k].second .resize(Nmfp, 0);
      }
      for (unsigned l=0; l<Nmfp; l++) {
	MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 50+Nphase+l, MPI_COMM_WORLD, &s);
	for (unsigned k=0; k<num; k++) tmp3[k].second[l] = tmp[k];
      }
      mfpI.insert(mfpI.end(), tmp3.begin(), tmp3.end());
    }
    
    // Sort the counters in prep for quantiles
    // 
    std::sort(ratI.begin(),  ratI.end()); 
    std::sort(denI.begin(),  denI.end()); 
    std::sort(volI.begin(),  volI.end()); 
    std::sort(temI.begin(),  temI.end()); 
    std::sort(delI.begin(),  delI.end()); 
    std::sort(selI.begin(),  selI.end()); 
    std::sort(phsI.begin(),  phsI.end());
    std::sort(mfpI.begin(),  mfpI.end());
    std::sort(kerI.begin(),  kerI.end());
    std::sort(derI.begin(),  derI.end());
    
    collnum = selI.size();
    coolnum = kerI.size();
    
    mfp_  .resize(quantiles.size());
    ts_   .resize(quantiles.size());
    coll_ .resize(quantiles.size());
    cool_ .resize(quantiles.size());
    rate_ .resize(quantiles.size());

    for (unsigned j=0; j<quantiles.size(); j++) {
      if (mfpI.size())
	mfp_[j]  = mfpI [Qi(quantiles[j],mfpI.size())].first;
      else
	mfp_[j] = 0;
      if (ratI.size())
	ts_[j]   = ratI [Qi(quantiles[j],ratI.size()) ];
      else
	ts_[j]   = 0;
      if (selI.size())
	coll_[j] = selI [Qi(quantiles[j],selI.size()) ];
      else
	coll_[j] = 0;
      if (kerI.size())
	cool_[j] = ratI [Qi(quantiles[j],kerI.size()) ];
      else
	cool_[j] = 0;
      if (derI.size())
	rate_[j] = derI [Qi(quantiles[j],derI.size()) ];
      else
	rate_[j] = 0;
    }
    
    if (SORTED) {
      ostringstream ostr;
      ostr << outdir << runtag << ".collide." << this_step;
      ofstream out(ostr.str().c_str());
      out << left << setw(8) << "# N" // Header
	  << setw(18) << "| MFP/L"
	  << setw(18) << "| Cyl radius (MFP)"
	  << setw(18) << "| Vertical (MFP)"
	  << setw(18) << "| Sph radius (MFP)"
	  << setw(18) << "| Density(MFP)"
	  << setw(18) << "| Volume(MFP)"
	  << setw(18) << "| TOF/TS"
	  << setw(18) << "| Density"
	  << setw(18) << "| Cell vol"
	  << setw(18) << "| Cell temp"
	  << setw(18) << "| Cool/part"
	  << setw(18) << "| Number/Nsel"
	  << endl;
      out << "# " << setw(6) << 1;
      for (unsigned k=2; k<13; k++) out << "| " << setw(16) << k;
      out << endl;
      cout << "SORTED: " << mfpI.size() << " cells" << endl;
      for (unsigned j=0; j<mfpI.size(); j++)
	out << setw(8) << j 
	    << setw(18) << mfpI[j].first
	    << setw(18) << mfpI[j].second[0]
	    << setw(18) << mfpI[j].second[1]
	    << setw(18) << mfpI[j].second[2]
	    << setw(18) << mfpI[j].second[3]
	    << setw(18) << mfpI[j].second[4]
	    << setw(18) << ratI[j] 
	    << setw(18) << denI[j] 
	    << setw(18) << volI[j] 
	    << setw(18) << temI[j] 
	    << setw(18) << delI[j] 
	    << setw(18) << selI[j] 
	    << endl;
      
      out << flush;
      out.close();
    }
    
    
    if (PHASE) {
      ostringstream ostr;
      ostr << outdir << runtag << ".phase." << this_step;
      ofstream out(ostr.str().c_str());
      out << left << setw(8) << "# N" // Header
	  << setw(18) << "| Density"
	  << setw(18) << "| Temp"
	  << setw(18) << "| Number"
	  << setw(18) << "| Mass"
	  << setw(18) << "| Volume"
	  << endl;
      out << "# " << setw(6) << 1;
      for (unsigned k=2; k<7; k++) out << "| " << setw(16) << k;
      out << endl;
      for (unsigned j=0; j<phsI.size(); j++) {
	out << setw(8) << j << setw(18) << phsI[j].first;
	for (unsigned k=0; k<Nphase; k++) 
	  out << setw(18) << phsI[j].second[k];
	out << endl;
      }
    }
    
  } else {
    unsigned num = tmfpst.size();
    unsigned nmt = tdelt.size();
    unsigned nmb = derat.size();
    
    MPI_Send(&nmt, 1, MPI_UNSIGNED, 0, 37, MPI_COMM_WORLD);
    MPI_Send(&nmb, 1, MPI_UNSIGNED, 0, 38, MPI_COMM_WORLD);
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&tsrat[0],  num, MPI_DOUBLE, 0, 40, MPI_COMM_WORLD);
    MPI_Send(&tdens[0],  num, MPI_DOUBLE, 0, 41, MPI_COMM_WORLD);
    MPI_Send(&tvolc[0],  num, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
    MPI_Send(&ttemp[0],  nmt, MPI_DOUBLE, 0, 43, MPI_COMM_WORLD);
    MPI_Send(&tdelt[0],  nmt, MPI_DOUBLE, 0, 44, MPI_COMM_WORLD);
    MPI_Send(&tseln[0],  num, MPI_DOUBLE, 0, 45, MPI_COMM_WORLD);
    MPI_Send(&kerat[0],  nmb, MPI_DOUBLE, 0, 46, MPI_COMM_WORLD);
    MPI_Send(&derat[0],  nmb, MPI_DOUBLE, 0, 47, MPI_COMM_WORLD);
    
    vector<double> tmt(nmt), tmp(num);
    
    for (unsigned k=0; k<nmt; k++) tmt[k] = tphase[k].first;
    MPI_Send(&tmt[0], nmt, MPI_DOUBLE, 0, 48, MPI_COMM_WORLD);
    for (unsigned l=0; l<Nphase; l++) {
      for (unsigned k=0; k<nmt; k++) tmt[k] = tphase[k].second[l];
      MPI_Send(&tmt[0], nmt, MPI_DOUBLE, 0, 49+l, MPI_COMM_WORLD);
    }
    
    for (unsigned k=0; k<num; k++) tmp[k] = tmfpst[k].first;
    MPI_Send(&tmp[0], num, MPI_DOUBLE, 0, 49+Nphase, MPI_COMM_WORLD);
    for (unsigned l=0; l<Nmfp; l++) {
      for (unsigned k=0; k<num; k++) tmp[k] = tmfpst[k].second[l];
      MPI_Send(&tmp[0], num, MPI_DOUBLE, 0, 50+Nphase+l, MPI_COMM_WORLD);
    }
  }
}

void Collide::EPSM(pCell* const cell, int id)
{
  if (cell->bods.size()<2) return;
  
#ifdef USE_GPTL
  GPTLstart("Collide::EPSM");
#endif
  // Compute mean and variance in each dimension
  // 
  EPSMT[id][0].start();
  
  vector<double> mvel(3, 0.0), disp(3, 0.0);
  double mass = 0.0;
  double Exes = 0.0;
  unsigned nbods = cell->bods.size();
  double coolheat = getCoolingRate(id);
  
  for (auto ib : cell->bods) {
    
    Particle* const p = tree->Body(ib);
    if (p->mass<=0.0 || std::isnan(p->mass)) {
      cout << "[crazy mass]";
    }
    for (unsigned k=0; k<3; k++) {
      mvel[k] += p->mass*p->vel[k];
      disp[k] += p->mass*p->vel[k]*p->vel[k];
      if (fabs(p->pos[k])>1.0e6 || std::isnan(p->pos[k])) {
	cout << "[crazy pos]";
      }
      if (fabs(p->vel[k])>1.0e6 || std::isnan(p->vel[k])) {
	cout << "[crazy vel]";
      }
    }
    mass += p->mass;
    
    // Compute the total  undercooled (-) or overcooled (+) energy.
    // That is Exes must be added to the internal energy before cooling
    // at this step.  If use_exes<0, Exes will remain zero (0).
    if (use_exes>=0) {
      Exes += p->dattrib[use_exes];
      p->dattrib[use_exes] = 0;
    }
  }
  EPSMTSoFar[id][0] = EPSMT[id][0].stop();
  
  //
  // Can't do anything if the gas has no mass
  //
  if (mass<=0.0) return;
  
  //
  // Compute the thermal (e.g. internal) energy
  //
  EPSMT[id][1].start();
  
  double Einternal = 0.0, Enew;
  for (unsigned k=0; k<3; k++) {
    mvel[k] /= mass;
    // Disp is variance here
    disp[k] = (disp[k] - mvel[k]*mvel[k]*mass)/mass;
    
    // Crazy value?
    if (disp[k]<0.0) disp[k] = 0.0;
    
    // Total kinetic energy in COV frame
    Einternal += 0.5*mass*disp[k];
  }
  
  EPSMTSoFar[id][1] = EPSMT[id][1].stop();
  
  //
  // Can't collide if with no internal energy
  //
  if (Einternal<=0.0) return;
  
  //
  // Correct 1d vel. disp. after cooling
  // 
  EPSMT[id][2].start();
  
  double Emin = 1.5*boltz*TFLOOR * mass/mp * 
    TreeDSMC::Munit/TreeDSMC::Eunit;
  
  // Einternal+Exes is the amount that
  // is *really* available for cooling
  // after excess or deficit is included
  //
  // Exes will be - if more energy still
  // needs to be removed and + if too much
  // energy was removed by cooling last step
  // 
  
  // Again, this should be moved to the
  // derived class
  
  if (Einternal + Exes - Emin > coolheat) {
    Enew = Einternal + Exes - coolheat;
  } else {
    Enew = min<double>(Emin, Einternal);
    
    decelT[id] += Einternal - Enew + Exes - coolheat;
    
    if (TSDIAG) {
      if (coolheat-Exes>0.0) {
	
	int indx = (int)floor(log(Einternal/coolheat) /
			      (log(2.0)*TSPOW) + 5);
	if (indx<0 ) indx = 0;
	if (indx>10) indx = 10;
	
	EoverT[id][indx] += mass;
      }
    }
  }
  // Compute the mean 1d vel.disp. from the
  // new internal energy value
  double mdisp = sqrt(Enew/mass/3.0);
  
  EPSMTSoFar[id][2] = EPSMT[id][2].stop();
  
  // Sanity check
  // 
  if (mdisp<=0.0 || std::isnan(mdisp) || std::isinf(mdisp)) {
    cout << "Process " << myid  << " id " << id 
	 << ": crazy values, mdisp=" << mdisp << " Enew=" << Enew
	 << " Eint=" << Einternal << " nbods=" << nbods << endl;
    return;
  }
  // Realize new velocities for all particles
  // 
  if (PULLIN) {
    EPSMT[id][3].start();
    
    double R=0.0, T=0.0;	// [Shuts up the compile-time warnings]
    const double sqrt3 = sqrt(3.0);
    
    if (nbods==2) {
      Particle* p1 = tree->Body(cell->bods[0]);
      Particle* p2 = tree->Body(cell->bods[1]);
      for (unsigned k=0; k<3; k++) {
	R = unit(random_gen);
	if (unit(random_gen)>0.5)
	  p1->vel[k] = mvel[k] + mdisp;
	else 
	  p1->vel[k] = mvel[k] - mdisp;
	p2->vel[k] = 2.0*mvel[k] - p1->vel[k];
      }
      
    } else if (nbods==3) {
      Particle* p1 = tree->Body(cell->bods[0]);
      Particle* p2 = tree->Body(cell->bods[1]);
      Particle* p3 = tree->Body(cell->bods[2]);
      double v2, v3;
      for (unsigned k=0; k<3; k++) {
	T = 2.0*M_PI*unit(random_gen);
	v2 = M_SQRT2*mdisp*cos(T);
	v3 = M_SQRT2*mdisp*sin(T);
	p1->vel[k] = mvel[k] - M_SQRT2*v2/sqrt3;
	p2->vel[k] = p1->vel[k] + (sqrt3*v2 - v3)/M_SQRT2;
	p3->vel[k] = p2->vel[k] + M_SQRT2*v3;
      }
    } else if (nbods==4) {
      Particle* p1 = tree->Body(cell->bods[0]);
      Particle* p2 = tree->Body(cell->bods[1]);
      Particle* p3 = tree->Body(cell->bods[2]);
      Particle* p4 = tree->Body(cell->bods[3]);
      double v2, v3, e2, e4, v4;
      for (unsigned k=0; k<3; k++) {
	R = unit(random_gen);
	e2 = mdisp*mdisp*(1.0 - R*R);
	T = 2.0*M_PI*unit(random_gen);
	v2 = sqrt(2.0*e2)*cos(T);
	v3 = sqrt(2.0*e2)*sin(T);
	p1->vel[k] = mvel[k] - sqrt3*v2/2.0;
	p2->vel[k] = p1->vel[k] + (2.0*v2 - M_SQRT2*v3)/sqrt3;
	e4 = mdisp*mdisp*R*R;
	if (unit(random_gen)>0.5) v4 =  sqrt(2.0*e4);
	else               v4 = -sqrt(2.0*e4);
	p3->vel[k] = p2->vel[k] + (sqrt3*v3 - v4)/M_SQRT2;
	p4->vel[k] = p3->vel[k] + M_SQRT2*v4;
      }
      
    } else {
      
      Particle *Pm1, *P00, *Pp1;
      vector<double> Tk, v(nbods), e(nbods);
      int kmax, dim, jj;
      bool Even = (nbods/2*2 == nbods);
      
      for (int k=0; k<3; k++) {

	if (Even) { 
	  // Even
	  kmax = nbods;
	  dim = kmax/2-1;
	  Tk .resize(dim);
	  for (int m=0; m<dim; m++) 
	    Tk[m] = pow(unit(random_gen), 1.0/(kmax/2 - m - 1.5));
	} else {			
	  // Odd
	  kmax = nbods-1;
	  dim = kmax/2-1;
	  Tk .resize(dim);
	  for (int m=0; m<dim; m++) 
	    Tk[m] = pow(unit(random_gen), 1.0/(kmax/2 - m - 1.0));
	}
	
	e[1] = mdisp*mdisp*(1.0 - Tk[0]);
	T = 2.0*M_PI*unit(random_gen);
	v[1] = sqrt(2.0*e[1])*cos(T);
	v[2] = sqrt(2.0*e[1])*sin(T);
	
	P00 = tree->Body(cell->bods[0]);
	Pp1 = tree->Body(cell->bods[1]);
	
	P00->vel[k] = mvel[k] - sqrt(nbods-1)*v[1]/sqrt(nbods);
	Pp1->vel[k] = P00->vel[k] + (sqrt(nbods)*v[1] - sqrt(nbods-2)*v[2])/sqrt(nbods-1);
	
	double prod = 1.0;
	for (int j=4; j<kmax-1; j+=2) {
	  jj = j-1;
	  
	  Pm1 = tree->Body(cell->bods[jj-2]);
	  P00 = tree->Body(cell->bods[jj-1]);
	  Pp1 = tree->Body(cell->bods[jj  ]);
	  
	  prod *= Tk[j/2-2];
	  e[jj] = mdisp*mdisp*(1.0 - Tk[j/2-1])*prod;
	  T = 2.0*M_PI*unit(random_gen);
	  v[jj]   = sqrt(2.0*e[jj])*cos(T);
	  v[jj+1] = sqrt(2.0*e[jj])*sin(T);
	  
	  P00->vel[k] = Pm1->vel[k] + 
	    (sqrt(3.0+nbods-j)*v[jj-1] - sqrt(1.0+nbods-j)*v[jj]  )/sqrt(2.0+nbods-j);
	  Pp1->vel[k] = P00->vel[k] +
	    (sqrt(2.0+nbods-j)*v[jj  ] - sqrt(    nbods-j)*v[jj+1])/sqrt(1.0+nbods-j);
	}
	
	prod *= Tk[kmax/2-2];
	e[kmax-1] = mdisp*mdisp*prod;
	
	if (Even) {
	  if (unit(random_gen)>0.5) v[nbods-1] =  sqrt(2.0*e[kmax-1]);
	  else               v[nbods-1] = -sqrt(2.0*e[kmax-1]);
	} else {
	  T = 2.0*M_PI*unit(random_gen);
	  v[nbods-2] = sqrt(2.0*e[kmax-1])*cos(T);
	  v[nbods-1] = sqrt(2.0*e[kmax-1])*sin(T);
	  
	  Pm1 = tree->Body(cell->bods[nbods-4]);
	  P00 = tree->Body(cell->bods[nbods-3]);
	  
	  P00->vel[k] = Pm1->vel[k] + (2.0*v[nbods-3] - M_SQRT2*v[nbods-2])/sqrt3;
	}
	
	Pm1 = tree->Body(cell->bods[nbods-3]);
	P00 = tree->Body(cell->bods[nbods-2]);
	Pp1 = tree->Body(cell->bods[nbods-1]);
	
	P00->vel[k] = Pm1->vel[k] + (sqrt3*v[nbods-2] - v[nbods-1])/M_SQRT2;
	Pp1->vel[k] = P00->vel[k] + M_SQRT2*v[nbods-1];
      }
    }
    
    // End Pullin algorithm
    
    EPSMTSoFar[id][1] = EPSMT[id][3].stop();
    
  } else {
    
    EPSMT[id][4].start();
    
    //
    // Realize a distribution with internal dispersion only
    //
    vector<double> Tmvel(3, 0.0);
    vector<double> Tdisp(3, 0.0);
    
    for (unsigned j=0; j<nbods; j++) {
      Particle* p = tree->Body(cell->bods[j]);
      
      for (unsigned k=0; k<3; k++) {
	p->vel[k] = mdisp*norm(random_gen);
	Tmvel[k] += p->mass*p->vel[k];
	Tdisp[k] += p->mass*p->vel[k]*p->vel[k];
	if (fabs(p->vel[k])>1e6 || std::isnan(p->vel[k])) {
	  cout << "[Collide crazy vel indx=" << p->indx 
	       << hex << ", key=" << p->key << dec << "]";
	}
      }
    }
    
    //
    // Compute mean and variance
    // 
    double Tmdisp = 0.0;
    for (unsigned k=0; k<3; k++) {
      Tmvel[k] /= mass;
      Tdisp[k] = (Tdisp[k] - Tmvel[k]*Tmvel[k]*mass)/mass;
      Tmdisp += Tdisp[k];
    }
    Tmdisp = sqrt(0.5*Tmdisp/3.0);
    
    //
    // Sanity check
    // 
    if (Tmdisp<=0.0 || std::isnan(Tmdisp) || std::isinf(Tmdisp)) {
      cout << "Process " << myid  << " id " << id 
	   << ": crazy values, Tmdisp=" << Tmdisp << " mdisp=" << mdisp 
	   << " nbods=" << nbods << endl;
      return;
    }
    
    //
    // Enforce energy and momentum conservation
    // 
    for (unsigned j=0; j<nbods; j++) {
      Particle* p = tree->Body(cell->bods[j]);
      for (unsigned k=0; k<3; k++)
	p->vel[k] = mvel[k] + (p->vel[k]-Tmvel[k])*mdisp/Tmdisp;
    }
    
    EPSMTSoFar[id][4] = EPSMT[id][4].stop();
  }
  
  EPSMT[id][5].start();
  
  //
  // Debugging sanity check
  // 
  if (0) {
    vector<double> mvel1(3, 0.0), disp1(3, 0.0);
    double dvel = 0.0;
    double mass1 = 0.0;
    double Efinal = 0.0;
    double mdisp1 = 0.0;
    
    for (unsigned j=0; j<nbods; j++) {
      Particle* p = tree->Body(cell->bods[j]);
      for (unsigned k=0; k<3; k++) {
	mvel1[k] += p->mass*p->vel[k];
	disp1[k] += p->mass*p->vel[k]*p->vel[k];
      }
      mass1 += p->mass;
    }
    
    for (unsigned k=0; k<3; k++) {
      mvel1[k] /= mass1;
      // Disp is variance here
      disp1[k] = (disp1[k] - mvel1[k]*mvel1[k]*mass1)/mass1;
      mdisp1 += disp1[k];
      
      // Crazy value?
      if (disp1[k]<0.0) disp1[k] = 0.0;
      
      // Total kinetic energy in COV frame
      Efinal += 0.5*mass1*disp1[k];
      dvel += (mvel[k] - mvel1[k])*(mvel[k] - mvel1[k]);
    }
    
    if (fabs(Efinal - Einternal)>1.0e-8*Einternal) {
      cerr << "Process " << myid << ": Collide::EPSM: energy boo-boo,"
	   << "  nbods="    << nbods
	   << "  Efinal="   << Efinal
	   << "  Einter="   << Einternal
	   << "  mass="     << mass1
	   << "  delta D="  << (mdisp - sqrt(mdisp1/3.0))/mdisp
	   << "  delta E="  << Enew - Einternal
	   << "  delta F="  << (Efinal - Einternal)/Einternal
	   << "  delta V="  << sqrt(dvel)
	   << ")"
	   << endl;
    }
  }
  
  //
  // Record diagnostics
  // 
  lostSoFar_EPSM[id] += Einternal - Enew;
  epsm1T[id] += nbods;
  Nepsm1T[id]++;
  
#ifdef USE_GPTL
  GPTLstop("Collide::EPSM");
#endif
  
  EPSMTSoFar[id][5] = EPSMT[id][5].stop();
  
}


void Collide::EPSMtimingGather()
{
  for (int i=0; i<nEPSMT; i++) EPSMtime[i] = 0;
  
  for (int n=0; n<nthrds; n++) {
    for (int i=0; i<nEPSMT; i++) {
      EPSMtime[i] += EPSMTSoFar[n][i];
      EPSMT[n][i].reset();
    }
  }
  
  if (myid==0) {
    MPI_Reduce(MPI_IN_PLACE, &EPSMtime[0], nEPSMT, MPI_LONG, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&EPSMtime[0], MPI_IN_PLACE, nEPSMT, MPI_LONG, MPI_SUM, 0,
	       MPI_COMM_WORLD);
  }
}

void Collide::EPSMtiming(ostream& out)
{
  if (EPSMratio<=0.0) return;
  
  if (myid==0) {
    
    const char labels[nEPSMT][20] = {"Mean/Var",
				     "Energy",
				     "1-d vel",
				     "Pullin",
				     "Standard",
				     "Final"};
    long sum = 0;
    for (int i=0; i<nEPSMT; i++) sum += EPSMtime[i];
    
    if (sum>0.0) {
      out << setfill('-') << setw(40) << '-' << setfill(' ') << endl
	  << "EPSM timing" << endl
	  << "-----------" << endl;
      for (int i=0; i<nEPSMT; i++)
	out << left << setw(20) << labels[i] << setw(10) << EPSMtime[i] 
	    << setw(10) << fixed << 100.0*EPSMtime[i]/sum << "%" << endl;
      out << setfill('-') << setw(40) << '-' << setfill(' ') << endl;
    }      
  }
}

void Collide::list_sizes()
{
  string sname = outdir + runtag + ".collide_storage";
  for (int n=0; n<numprocs; n++) {
    if (myid==n) {
      ofstream out(sname.c_str(), ios::app);
      if (out) {
	out << setw(18) << tnow
	    << setw(6)  << myid;
	list_sizes_proc(&out);
	out << endl;
	if (myid==numprocs-1) out << endl;
	out.close();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void Collide::list_sizes_proc(ostream* out)
{
  *out << setw(12) << numcnt.size()
       << setw(12) << colcnt.size()
       << setw(12) << tsrat.size()
       << setw(12) << derat.size()
       << setw(12) << tdens.size()
       << setw(12) << tvolc.size()
       << setw(12) << ttemp.size()
       << setw(12) << tdelt.size()
       << setw(12) << tseln.size()
       << setw(12) << tphase.size()
       << setw(12) << (tphaseT.size() ? tphaseT[0].size() : (size_t)0)
       << setw(12) << tmfpst.size()
       << setw(12) << (tmfpstT.size() ? tmfpstT[0].size() : (size_t)0)
       << setw(12) << (numcntT.size() ? numcntT[0].size() : (size_t)0)
       << setw(12) << (colcntT.size() ? colcntT[0].size() : (size_t)0)
       << setw(12) << error1T.size()
       << setw(12) << col1T.size()
       << setw(12) << epsm1T.size()
       << setw(12) << Nepsm1T.size()
       << setw(12) << KEtotT.size()
       << setw(12) << KElostT.size()
       << setw(12) << tmassT.size()
       << setw(12) << decelT.size()
       << setw(12) << (mfpratT.size() ? mfpratT[0].size() : (size_t)0)
       << setw(12) << (tsratT.size() ? tsratT[0].size() : (size_t)0)
       << setw(12) << (tdensT.size() ? tdensT[0].size() : (size_t)0)
       << setw(12) << (tvolcT.size() ? tvolcT[0].size() : (size_t)0)
       << setw(12) << (ttempT.size() ? ttempT[0].size() : (size_t)0)
       << setw(12) << (tselnT.size() ? tselnT[0].size() : (size_t)0)
       << setw(12) << (deratT.size() ? deratT[0].size() : (size_t)0)
       << setw(12) << (tdeltT.size() ? tdeltT[0].size() : (size_t)0)
       << setw(12) << (tdispT.size() ? tdispT[0].size() : (size_t)0)
       << setw(12) << tdiag.size()
       << setw(12) << tdiag1.size()
       << setw(12) << tdiag0.size()
       << setw(12) << tcool.size()
       << setw(12) << tcool1.size()
       << setw(12) << tcool0.size()
       << setw(12) << (cellist.size() ? cellist[0].size() : (size_t)0)
       << setw(12) << disptot.size()
       << setw(12) << lostSoFar_EPSM.size();
}


void Collide::CollectTiming()
{
  int nf = 5, c;
  if (TIMING) nf += 11;
  
  vector<double> in(nf, 0.0);
  vector< vector<double> > out(3);
  for (int i=0; i<3; i++) out[i] .resize(nf);
  
  in[0] += forkSoFar;
  in[1] += snglSoFar;
  in[2] += waitSoFar;
  in[3] += diagSoFar;
  in[4] += joinSoFar;
  
  for (int n=0; n<nthrds; n++) {
    c = 5;
    in[c++] += listSoFar[n];
    in[c++] += initSoFar[n];
    in[c++] += collSoFar[n];
    in[c++] += elasSoFar[n];
    in[c++] += cellSoFar[n];
    in[c++] += epsmSoFar[n];
    in[c++] += coolSoFar[n];
    in[c++] += stat1SoFar[n];
    in[c++] += stat2SoFar[n];
    in[c++] += stat3SoFar[n];
    in[c++] += collCnt[n];
  }
  
  // Minimum
  MPI_Reduce(&in[0], &out[0][0], nf, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  
  // Sum
  MPI_Reduce(&in[0], &out[1][0], nf, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // Maximum
  MPI_Reduce(&in[0], &out[2][0], nf, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  
  c = 0;
  forkSum[0]  = out[0][c];
  forkSum[1]  = out[1][c]/numprocs;
  forkSum[2]  = out[2][c];
  
  c++;
  snglSum[0]  = out[0][c];
  snglSum[1]  = out[1][c]/numprocs;
  snglSum[2]  = out[2][c];
  
  c++;
  waitSum[0]  = out[0][c];
  waitSum[1]  = out[1][c]/numprocs;
  waitSum[2]  = out[2][c];
  
  c++;
  diagSum[0]  = out[0][c];
  diagSum[1]  = out[1][c]/numprocs;
  diagSum[2]  = out[2][c];
  
  c++;
  joinSum[0]  = out[0][c];
  joinSum[1]  = out[1][c]/numprocs;
  joinSum[2]  = out[2][c];
  
  if (TIMING) {
    
    c++;
    listSum[0]  = out[0][c];
    listSum[1]  = out[1][c]/numprocs;
    listSum[2]  = out[2][c];
    
    c++;
    initSum[0]  = out[0][c];
    initSum[1]  = out[1][c]/numprocs;
    initSum[2]  = out[2][c];
    
    c++;
    collSum[0]  = out[0][c];
    collSum[1]  = out[1][c]/numprocs;
    collSum[2]  = out[2][c];
    
    c++;
    elasSum[0]  = out[0][c];
    elasSum[1]  = out[1][c]/numprocs;
    elasSum[2]  = out[2][c];
    
    c++;
    cellSum[0]  = out[0][c];
    cellSum[1]  = out[1][c]/numprocs;
    cellSum[2]  = out[2][c];
    
    c++;
    epsmSum[0]  = out[0][c];
    epsmSum[1]  = out[1][c]/numprocs;
    epsmSum[2]  = out[2][c];
    
    c++;
    coolSum[0]  = out[0][c];
    coolSum[1]  = out[1][c]/numprocs;
    coolSum[2]  = out[2][c];
    
    c++;
    stat1Sum[0] = out[0][c];
    stat1Sum[1] = out[1][c]/numprocs;
    stat1Sum[2] = out[2][c];
    
    c++;
    stat2Sum[0] = out[0][c];
    stat2Sum[1] = out[1][c]/numprocs;
    stat2Sum[2] = out[2][c];
    
    c++;
    stat3Sum[0] = out[0][c];
    stat3Sum[1] = out[1][c]/numprocs;
    stat3Sum[2] = out[2][c];
    
    c++;
    numbSum[0]  = out[0][c];
    numbSum[1]  = floor(out[1][c]/numprocs+0.5);
    numbSum[2]  = out[2][c];
  }
  
  // Reset the timers
  
  forkTime.reset();
  snglTime.reset();
  waitTime.reset();
  joinTime.reset();
  diagTime.reset();
  
  for (int n=0; n<nthrds; n++) {
    listTime [n].reset(); 
    initTime [n].reset(); 
    collTime [n].reset(); 
    stat1Time[n].reset();
    stat2Time[n].reset();
    stat3Time[n].reset();
    coolTime [n].reset(); 
    elasTime [n].reset(); 
    cellTime [n].reset(); 
    epsmTime [n].reset(); 
    collCnt  [n] = 0;
  }
  
}


template<typename T>
void Collide::colldeHelper(ostream& out, const char* lab, vector<T>& v)
{
  out << left << fixed << setw(18) << lab 
      << right << "[" << setprecision(6)
      << setw(10) << v[0] << ", " 
      << setw(10) << v[1] << ", "
      << setw(10) << v[2] << "]" << endl << left;
}

void Collide::colldeTime(ostream& out) 
{ 
  out << "-----------------------------------------------------" << endl;
  out << "-----Collide timing----------------------------------" << endl;
  out << "-----------------------------------------------------" << endl;
  
  colldeHelper<double>(out, "Thread time",  forkSum); 
  colldeHelper<double>(out, "Joined time",  snglSum);
  colldeHelper<double>(out, "Waiting time", waitSum);
  colldeHelper<double>(out, "Join time",    joinSum);
  colldeHelper<double>(out, "Diag time",    diagSum);
  
  out << left 
      << setw(18) << "Body count"   << bodycount    
      << " [" << bodycount/stepcount << "]" << endl
      << setw(18) << "Step count"   << stepcount    << endl
      << endl;
  
  if (TIMING) {
    colldeHelper<double>(out, "All cells",  cellSum); 
    colldeHelper<double>(out, "List bods",  listSum);
    colldeHelper<double>(out, "Init",       initSum);
    colldeHelper<double>(out, "Collide",    collSum); 
    colldeHelper<double>(out, "Inelastic",  elasSum); 
    colldeHelper<double>(out, "EPSM",       epsmSum); 
    colldeHelper<double>(out, "Stat#1",     stat1Sum); 
    colldeHelper<double>(out, "Stat#2",     stat2Sum); 
    colldeHelper<double>(out, "Stat#3",     stat3Sum); 
    colldeHelper<int   >(out, "Cell count", numbSum); 
    out << endl;
  }
  
  stepcount = 0;
  bodycount = 0;
  
  out << "-----------------------------------------------------" << endl;
}

void Collide::tsdiag(ostream& out) 
{
  if (!TSDIAG) return;
  
  if (tdiag.size()==numdiag) {
    
    out << "-----------------------------------------------------" << endl;
    out << "-----Time step diagnostics---------------------------" << endl;
    out << "-----------------------------------------------------" << endl;
    out << right << setw(8) << "2^n" << setw(15) << "TS ratio"
	<< setw(15) << "Size/Vel";
    if (use_delt>=0) out << setw(15) << "Kinetic/Cool";
    out << endl << setprecision(3);
    out << "-----------------------------------------------------" << endl;
    for (unsigned k=0; k<numdiag; k++) {
      double rat = pow(4.0, -5.0+k);
      out << setw(8)  << -10+2*static_cast<int>(k)
	  << (((rat<1.0e-02 || rat>1.0e+06) && rat>0.0) ? scientific : fixed)
	  << setw(15) << rat
	  << setw(15) << tdiag[k];
      if (use_delt>=0) out << setw(15) << tcool[k];
      out << endl;
      tdiag[k] = 0;
      if (use_delt>=0) tcool[k] = 0;
    }
    out << "-----------------------------------------------------" << endl;
    out << left;
  }
  
  
  if (Eover.size()==numdiag) {
    
    double emass = 0.0;
    for (unsigned k=0; k<numdiag; k++) emass += Eover[k];
    if (emass>0.0) {
      out << "-----Cooling rate diagnostics------------------------" << endl;
      out << "-----------------------------------------------------" << endl;
      out << right << setw(8) << "2^n" << setw(15) << "Ratio"
	  << setw(15) << "KE/Cool(%)" << endl;
      out << "-----------------------------------------------------" << endl;
      
      for (unsigned k=0; k<numdiag; k++) {
	double rat = pow(pow(2.0, TSPOW), -5.0+k);
	double val = Eover[k]*100.0/emass;
	out << setw(8)  << TSPOW*(-5 + static_cast<int>(k))
	    << (((rat<1.0e-02 || rat>1.0e+06) && rat>0.0) ? scientific : fixed)
	    << setw(15) << rat
	    << (((val<1.0e-02 || val>1.0e+06) && val>0.0) ? scientific : fixed)
	    << setw(15) << val << endl;
	Eover[k] = 0;
      }
      out << "-----------------------------------------------------" << endl;
    }
    out << left;
  }
  
  if (Cover.size()==numdiag) {
    
    double cmass = 0.0;
    for (unsigned k=0; k<numdiag; k++) cmass += Cover[k];
    if (cmass>0.0) {
      out << "-----CBA scale diagnostics--------------------------" << endl;
      out << "-----------------------------------------------------" << endl;
      out << right << setw(8) << "2^n" << setw(15) << "Ratio"
	  << setw(15) << "Diam/Side(%)" << endl;
      out << "-----------------------------------------------------" << endl;
      
      for (unsigned k=0; k<numdiag; k++) {
	double rat = pow(4.0, -5.0+k);
	double val = Cover[k]*100.0/cmass;
	out << setw(8)  << -10+2*static_cast<int>(k)
	    << (((rat<1.0e-02 || rat>1.0e+06) && rat>0.0) ? scientific : fixed)
	    << setw(15) << rat
	    << (((val<1.0e-02 || val>1.0e+06) && val>0.0) ? scientific : fixed)
	    << setw(15) << val << endl;
	Cover[k] = 0;
      }
      out << "-----------------------------------------------------" << endl;
    }
    out << left;
  }
  
}


void Collide::voldiag(ostream& out) 
{
  if (!VOLDIAG) return;
  
  if (Vcnt.size()==nbits) {
    // Find the smallest cell size
    unsigned nlast;
    for (nlast=nbits; nlast>0; nlast--)
      if (Vcnt[nlast-1]) break;
    
    out << "-----------------------------------------------------"
	<< "-----------------------------------------------------" << endl;
    out << "-----Volume cell diagnostics-------------------------"
	<< "-----------------------------------------------------" << endl;
    out << "-----------------------------------------------------"
      	<< "-----------------------------------------------------" << endl;
    out << right << setw(8) << "n" << setw(10) << "#"
	<< setw(12) << "Factor"
	<< setw(12) << "Density"   << setw(12) << "MFP/L"
	<< setw(12) << "Coll prob" << setw(12) << "Flight/L"
	<< setw(12) << "Particles" << setw(12) << "Root var";
    out << endl << setprecision(3) << scientific;
    out << "-----------------------------------------------------"
      	<< "-----------------------------------------------------" << endl;
    
    for (unsigned k=0; k<nlast; k++) {
      unsigned n  = k*nvold;
      unsigned nn = n + nvold - 1;
      double rat  = pow(2.0, -3.0*k);
      double nrm  = Vcnt[k] ? 1.0/Vcnt[k] : 0.0;
      out << setw(8) << k << setw(10) << Vcnt[k] << setw(12) << rat;
      for (unsigned l=0; l<nvold-1; l++) out << setw(12) << Vdbl[n+l]*nrm;
      // Variance term
      Vdbl[nn] = fabs(Vdbl[nn] - Vdbl[nn-1]*Vdbl[nn-1]*nrm);
      out << setw(12) << sqrt(Vdbl[nn]*nrm) << endl;
    }
    out << "-----------------------------------------------------"
	<< "-----------------------------------------------------" << endl;
    out << left;
  }
  
}


// For timestep computation

extern "C"
void *
tstep_thread_call(void *atp)
{
  thrd_pass_tstep *tp = (thrd_pass_tstep *)atp;
  Collide *p = (Collide *)tp->p;

  p->timestep_thread((void*)&tp->arg);

  return NULL;
}

void Collide::compute_timestep(double coolfrac)
{
  int errcode;
  void *retval;
  
  if (nthrds==1) {
    thrd_pass_tstep td;
    
    td.p            = this;
    td.arg.coolfrac = coolfrac;
    td.arg.id       = 0;
    
    tstep_thread_call(&td);
    
    return;
  }
  
  tdT = new thrd_pass_tstep [nthrds];
  t = new pthread_t [nthrds];
  
  if (!tdT) {
    cerr << "Process " << myid 
         << ": Collide::tstep_thread_call: error allocating memory for thread counters\n";
    exit(18);
  }
  if (!t) {
    cerr << "Process " << myid
         << ": Collide::tstep_thread_call: error allocating memory for thread\n";
    exit(18);
  }
  
  // Make the <nthrds> threads
  //
  for (int i=0; i<nthrds; i++) {
    tdT[i].p            = this;
    tdT[i].arg.coolfrac = coolfrac;
    tdT[i].arg.id       = i;
    
    errcode =  pthread_create(&t[i], 0, tstep_thread_call, &tdT[i]);
    if (errcode) {
      cerr << "Process " << myid;
      cerr << " Collide: cannot make thread " << i
	   << ", errcode=" << errcode << endl;
      exit(19);
    }
  }
  
  // Collapse the threads
  for (int i=0; i<nthrds; i++) {
    if ((errcode=pthread_join(t[i], &retval))) {
      cerr << "Process " << myid;
      cerr << " Collide::tstep_thread_call: thread join " << i
           << " failed, errcode=" << errcode << endl;
      exit(20);
    }
  }
  
  delete [] tdT;
  delete [] t;
  
  caller->print_timings("Collide: timestep thread timings", timer_list);
}


void * Collide::timestep_thread(void * arg)
{
  double coolfrac = (double)((tstep_pass_arguments*)arg)->coolfrac;
  int id          = (int)((tstep_pass_arguments*)arg)->id;
  
  thread_timing_beg(id);
  
  // Loop over cells, cell time-of-flight time for each particle
  //
  pCell *c;
  Particle *p;
  double L, DT, mscale;
  
  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {
    
    // Number of particles in this cell
    //
    c = cellist[id][j];
    L = c->Scale();
    
    // Look for collCount in cell attributes
    //
    std::map<std::string, int>::iterator itc = c->iattrib.find("collCount");
    double DTcoll = DBL_MAX;
    if (itc != c->iattrib.end()) {
				// Timestep for this cell
    double tau = dtime / (1<<c->maxplev);
				// If there have been collisions, set
				// time step based on scaling to
      if (itc->second > 0)	// target
	DTcoll = tau * collTnum / itc->second;
    }

    for (auto i : c->bods) {

      // Current particle
      //
      p = tree->Body(i);

      // Compute time of flight criterion
      //
      DT     = 1.0e40;
      mscale = 1.0e40;
      for (unsigned k=0; k<3; k++)
	mscale = std::min<double>(pHOT::sides[k]*L, mscale);
      
      // Size scale for multistep timestep calc.
      //
      p->scale = mscale;

      // Collision target
      //
      DT = std::min<double>(DT, DTcoll);

      // Compute cooling criterion timestep
      //
      if (use_delt>=0) {
	double v = p->dattrib[use_delt];
	if (v>0.0) DT = min<double>(DT, coolfrac*v);
      }

      p->dtreq = DT;
    }
  }
  
  thread_timing_end(id);
  
  return (NULL);
}

void Collide::energyExcess(double& ExesColl, double& ExesEPSM)
{
  // Sum up from all the threads
  // 
  for (int n=1; n<nthrds; n++) {
    exesCT[0] += exesCT[n];
    exesET[0] += exesET[n];
  }
  // Sum reduce result to root node
  // 
  MPI_Reduce(&exesCT[0], &ExesColl, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&exesET[0], &ExesEPSM, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // Zero out the thread accumulators
  // 
  for (int n=0; n<nthrds; n++) exesCT[0] = exesET[n] = 0.0;
}


void Collide::pre_collide_diag()
{
  // Clean thread variables
  diagTime.start();
  for (int n=0; n<nthrds; n++) {
    error1T[n] = 0;
    sel1T[n]   = 0;
    col1T[n]   = 0;
    epsm1T[n]  = 0;
    Nepsm1T[n] = 0;
    tmassT[n]  = 0;
    decolT[n]  = 0;
    decelT[n]  = 0;
    
    // For computing cell occupation #
    colcntT[n].clear();	// and collision counts
    numcntT[n].clear();
    
    // For computing MFP to cell size ratio 
    // and drift ratio
    if (MFPDIAG) {
      tsratT[n].clear();
      keratT[n].clear();
      deratT[n].clear();
      tdensT[n].clear();
      tvolcT[n].clear();
      ttempT[n].clear();
      tdeltT[n].clear();
      tselnT[n].clear();
      tphaseT[n].clear();
      tmfpstT[n].clear();
    }
    
    if (VOLDIAG) {
      for (unsigned k=0; k<nbits; k++) {
	VcntT[n][k] = 0;
	for (unsigned l=0; l<nvold; l++) VdblT[n][k*nvold+l] = 0.0;
      }
    }
    
    for (unsigned k=0; k<numdiag; k++) {
      if (TSDIAG) {
	tdiagT[n][k] = 0;
	EoverT[n][k] = 0;
      }
      if (use_delt>=0) tcoolT[n][k] = 0;
    }
    
    for (unsigned k=0; k<3; k++) tdispT[n][k] = 0;
  }
  if (TSDIAG) {
    for (unsigned k=0; k<numdiag; k++) tdiag1[k] = tdiag0[k] = 0;
    for (unsigned k=0; k<numdiag; k++) Eover1[k] = Eover0[k] = 0;
  }
  if (VOLDIAG) {
    for (unsigned k=0; k<nbits; k++) {
      Vcnt1[k] = Vcnt0[k] = 0;
      for (unsigned l=0; l<nvold; l++) 
	Vdbl1[k*nvold+l] = Vdbl0[k*nvold+l] = 0.0;
    }
  }
  
  if (use_delt>=0) 
    for (unsigned k=0; k<numdiag; k++) tcool1[k] = tcool0[k] = 0;
  
  diagTime.stop();
  
  if (CDEBUG) list_sizes();
}


unsigned Collide::post_collide_diag()
{
  unsigned sel=0, col=0;
  

  diagTime.start();
  // Diagnostics

  unsigned error1=0, error=0;
  
  unsigned sel1=0, col1=0;	// Count number of collisions
  unsigned epsm1=0, epsm=0, Nepsm1=0, Nepsm=0;
  
  // Dispersion test
  double mass1 = 0, mass0 = 0;
  vector<double> disp1(3, 0), disp0(3, 0);
  
  numcnt.clear();
  colcnt.clear();
  
  for (int n=0; n<nthrds; n++) {
    error1 += error1T[n];
    sel1   += sel1T[n];
    col1   += col1T[n];
    epsm1  += epsm1T[n];
    Nepsm1 += Nepsm1T[n];
    numcnt.insert(numcnt.end(), numcntT[n].begin(), numcntT[n].end());
    colcnt.insert(colcnt.end(), colcntT[n].begin(), colcntT[n].end());
    if (TSDIAG) {
      for (unsigned k=0; k<numdiag; k++) tdiag1[k] += tdiagT[n][k];
      for (unsigned k=0; k<numdiag; k++) Eover1[k] += EoverT[n][k];
    }
    if (VOLDIAG) {
      for (unsigned k=0; k<nbits; k++) {
	Vcnt1[k] += VcntT[n][k];
	for (unsigned l=0; l<nvold; l++)
	  Vdbl1[k*nvold+l] += VdblT[n][k*nvold+l];
      }
    }
    if (use_delt>=0) 
      for (unsigned k=0; k<numdiag; k++) tcool1[k] += tcoolT[n][k];
  }
  
  // For computing MFP to cell size ratio 
  // and drift ratio (diagnostic only)
  if (MFPDIAG) {
    tsrat.clear();
    kerat.clear();
    derat.clear();
    tdens.clear();
    tvolc.clear();
    ttemp.clear();
    tdelt.clear();
    tseln.clear();
    tphase.clear();
    tmfpst.clear();
    
    for (int n=0; n<nthrds; n++) {
      tsrat. insert(tsrat.end(),   tsratT[n].begin(),  tsratT[n].end());
      kerat. insert(kerat.end(),   keratT[n].begin(),  keratT[n].end());
      derat. insert(derat.end(),   deratT[n].begin(),  deratT[n].end());
      tdens. insert(tdens.end(),   tdensT[n].begin(),  tdensT[n].end());
      tvolc. insert(tvolc.end(),   tvolcT[n].begin(),  tvolcT[n].end());
      ttemp. insert(ttemp.end(),   ttempT[n].begin(),  ttempT[n].end());
      tdelt. insert(tdelt.end(),   tdeltT[n].begin(),  tdeltT[n].end());
      tseln. insert(tseln.end(),   tselnT[n].begin(),  tselnT[n].end());
      tphase.insert(tphase.end(), tphaseT[n].begin(), tphaseT[n].end());
      tmfpst.insert(tmfpst.end(), tmfpstT[n].begin(), tmfpstT[n].end());
    }
  }
  
  for (int n=0; n<nthrds; n++) {
    for (unsigned k=0; k<3; k++) disp1[k] += tdispT[n][k];
    mass1 += tmassT[n];
  }
  
  MPI_Reduce(&sel1, &sel, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&col1, &col, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&epsm1, &epsm, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Nepsm1, &Nepsm, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&error1, &error, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ncells, &numtot, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if (TSDIAG) {
    MPI_Reduce(&tdiag1[0], &tdiag0[0], numdiag, MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
    MPI_Reduce(&Eover1[0], &Eover0[0], numdiag, MPI_DOUBLE,   MPI_SUM, 0, 
	       MPI_COMM_WORLD);
  }
  if (VOLDIAG) {
    MPI_Reduce(&Vcnt1[0], &Vcnt0[0], nbits, MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
    MPI_Reduce(&Vdbl1[0], &Vdbl0[0], nbits*nvold, MPI_DOUBLE, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
  }
  if (use_delt>=0)
    MPI_Reduce(&tcool1[0], &tcool0[0], numdiag, MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
  MPI_Reduce(&disp1[0], &disp0[0], 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass1, &mass0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  seltot    += sel;
  coltot    += col;
  epsmtot   += epsm;
  epsmcells += Nepsm;
  errtot    += error;
  if (TSDIAG) {
    for (unsigned k=0; k<numdiag; k++) {
      tdiag[k] += tdiag0[k];
      Eover[k] += Eover0[k];
    }
  }
  if (VOLDIAG) {
    for (unsigned k=0; k<nbits; k++) {
      Vcnt[k] += Vcnt0[k];
      for (unsigned l=0; l<nvold; l++)
	Vdbl[k*nvold+l] += Vdbl0[k*nvold+l];
    }
  }
  if (use_delt>=0)
    for (unsigned k=0; k<numdiag; k++) tcool[k] += tcool0[k];
  for (unsigned k=0; k<3; k++) disptot[k] += disp0[k];
  masstot   += mass0;
  diagSoFar = diagTime.stop();
  
  return col;
}

void Collide::CPUHogGather()
{
  for (int i=0; i<2; i++) {
    CPUH[i]   = std::numeric_limits<long>::max();
    CPUH[6+i] = 0;
    CPUH[2+i] = CPUH[4+i] = CPUH[8+i] = CPUH[10+i] = -1;
  }
  
  for (int n=0; n<nthrds; n++) {
    for (int i=0; i<2; i++) {
      if (CPUH[i] > minUsage[n*2+i]) {
	CPUH[i]    = minUsage[n*2+i];
	CPUH[2+i]  = minPart [n*2+i];
	CPUH[4+i]  = minCollP[n*2+i];
      }
      if (CPUH[6+i] < maxUsage[n*2+i]) {
	CPUH[6+i]  = maxUsage[n*2+i];
	CPUH[8+i]  = maxPart [n*2+i];
	CPUH[10+i] = maxCollP[n*2+i];
      }
      // Clear values for next call
      minUsage[n*2+i] = std::numeric_limits<long>::max();
      maxUsage[n*2+i] = 0;
      minPart [n*2+i] = -1;
      maxPart [n*2+i] = -1;
      minCollP[n*2+i] = -1;
      maxCollP[n*2+i] = -1;
    }
  }
  
  vector<long> U(12*numprocs);
  
  MPI_Gather(&CPUH[0], 12, MPI_LONG, &U[0], 12, MPI_LONG, 0, MPI_COMM_WORLD);
  
  if (myid==0) {
    
    for (int i=0; i<2; i++) {
      CPUH[i]   = std::numeric_limits<long>::max();
      CPUH[6+i] = 0;
      CPUH[2+i] = CPUH[4+i] = CPUH[8+i] = CPUH[10+i] = -1;
    }
    
    for (int n=0; n<numprocs; n++) {
      for (int i=0; i<2; i++) {
	if (CPUH[i] > U[n*12+i]) {
	  CPUH[i]   = U[n*12+i];
	  CPUH[2+i] = U[n*12+2+i];
	  CPUH[4+i] = U[n*12+4+i];
	}
      }
      for (int i=6; i<8; i++) {
	if (CPUH[i] < U[n*12+i]) {
	  CPUH[i]   = U[n*12+i];
	  CPUH[2+i] = U[n*12+2+i];
	  CPUH[4+i] = U[n*12+4+i];
	}
      }
    }
  }
  
}

void Collide::CPUHog(ostream& out)
{
  if (myid==0) {
    const unsigned f = 8;
    out << "Extremal cell timing" << endl
	<< "--------------------" << endl
	<< "T=" << tnow << ",  mstep=" << mstep << endl << right;
    if (EPSMratio>0) {
      out << "                 "
	  << setw(f) << "DSMC"   << "  " << setw(f) << "EPSM"
	  << "      " 
	  << setw(f) << "DSMC"   << "  " << setw(f) << "EPSM"
	  << "      " 
	  << setw(f) << "DSMC"   << ", " << setw(f) << "EPSM" << endl
	  << "  Minimum usage=(" 
	  << setw(f) << CPUH[ 0] << ", " << setw(f) << CPUH[ 1]
	  << ")  N=(" 
	  << setw(f) << CPUH[ 2] << ", " << setw(f) << CPUH[ 3]
	  << ")  C=(" 
	  << setw(f) << CPUH[ 4] << ", " << setw(f) << CPUH[ 5] << ")" << endl
	  << "  Maximum usage=(" 
	  << setw(f) << CPUH[ 6] << ", " << setw(f) << CPUH[ 7]
	  << ")  N=("
	  << setw(f) << CPUH[ 8] << ", " << setw(f) << CPUH[ 9]
	  << ")  C=(" 
	  << setw(f) << CPUH[10] << ", " << setw(f) << CPUH[11] << ")" << endl;
    } else {
      out << endl
	  << "  Minimum usage=(" 
	  << setw(f) << CPUH[ 0] << ")  N=(" 
	  << setw(f) << CPUH[ 2] << ")  C=(" 
	  << setw(f) << CPUH[ 4] << ")" << endl
	  << "  Maximum usage=(" 
	  << setw(f) << CPUH[ 6] << ")  N=("
	  << setw(f) << CPUH[ 8] << ")  C=(" 
	  << setw(f) << CPUH[10] << ")" << endl;
    }
    out << endl;
  }
}

double Collide::hsDiameter()
{
  const double Bohr = 5.2917721092e-09;
  return hsdiam*Bohr*sqrt(crossfac)/TreeDSMC::Lunit;
}

void Collide::printSpecies(std::map<speciesKey, unsigned long>& spec,
			   double temp)
{
  if (myid) return;

  typedef std::map<speciesKey, unsigned long> spCountMap;
  typedef spCountMap::iterator spCountMapItr;

				// Field width
  const unsigned short wid = 16;

  std::ofstream dout;

				// Generate the file name
  if (species_file_debug.size()==0) {
    std::ostringstream sout;
    sout << outdir << runtag << ".species";
    species_file_debug = sout.str();

    // Check for existence of file
    //
    std::ifstream in (species_file_debug.c_str());
    
    // Write a new file?
    //
    if (in.fail()) {
      
      // Open the file for the first time
      //
      dout.open(species_file_debug.c_str());

      // Print the header
      //
      dout << "# " 
	   << std::setw(wid) << std::right << "Time "
	   << std::setw(wid) << std::right << "Temp ";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++) {
	std::ostringstream sout;
	sout << "(" << it->first.first << "," << it->first.second << ") ";
	dout << setw(wid) << right << sout.str();
      }
      dout << std::endl;
      
      dout << "# " 
	   << std::setw(wid) << std::right << "--------"
	   << std::setw(wid) << std::right << "--------";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
	dout << setw(wid) << std::right << "--------";
      dout << std::endl;
      
    }
  }

  // Open for append
  //
  if (!dout.is_open())
    dout.open(species_file_debug.c_str(), ios::out | ios::app);


  double tmass = 0.0;
  for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
    tmass += atomic_weights[it->first.first] * it->second;

				// Use total mass to print mass
				// fraction
  dout << "  " 
       << std::setw(wid) << std::right << tnow
       << std::setw(wid) << std::right << temp;

  for (spCountMapItr it=spec.begin(); it != spec.end(); it++) {
    if (tmass > 0.0) 
      dout << std::setw(wid) << std::right 
	   << atomic_weights[it->first.first] * it->second / tmass;
    else
      dout << std::setw(wid) << std::right << 0.0;
  }
  dout << std::endl;
}

void Collide::velocityUpdate(Particle* p1, Particle* p2, double cr)
{
  vector<double> vcm(3), vrel(3);

  // Center of mass velocity
  //
  double tmass = p1->mass + p2->mass;
  for(unsigned k=0; k<3; k++)
    vcm[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k]) / tmass;
	    
  double cos_th = 1.0 - 2.0*unit(random_gen);       // Cosine and sine of
  double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
  double phi    = 2.0*M_PI*unit(random_gen);	       // Collision angle phi
  
  vrel[0] = cr*cos_th;	  // Compute post-collision
  vrel[1] = cr*sin_th*cos(phi); // relative velocity
  vrel[2] = cr*sin_th*sin(phi);
  
  // Update post-collision velocities
  // 
  for(unsigned k=0; k<3; k++ ) {
    p1->vel[k] = vcm[k] + p2->mass/tmass*vrel[k];
    p2->vel[k] = vcm[k] - p1->mass/tmass*vrel[k];
  }
}


// Compute the mean molecular weight in atomic mass units
double Collide::molWeight()
{
  const double f_H = 0.76;	// Assume only hydrogen and helium here

  // Only do this once and cache the result
  if (mol_weight<0)  {
    mol_weight = 1.0 / (f_H/atomic_weights[1] + (1.0-f_H)/atomic_weights[2]);
  }

  return mol_weight;
}

void Collide::NTCgather()
{
  if (NTC and DEBUG_NTC) {

    qq.clear();

    unsigned Ovr=0, Acc=0, Tot=0, Max=0, gbMax;
    for (int n=0; n<nthrds; n++) {
      Ovr += ntcOvr[n];
      Acc += ntcAcc[n];
      Tot += ntcTot[n];
      Max  = std::max<unsigned>(ntcAcc[n], Max);
				// Reset the counters
      ntcOvr[n] = ntcAcc[n] = ntcTot[n] = 0;
    }

				// Accumulate into id=0
    ntcSum.clear();
    for (int n=0; n<nthrds; n++) {
      ntcSum.insert(ntcSum.end(), ntcVal[n]->begin(), ntcVal[n]->end());
    }

				// Accumulate into wgtTot
    wgtSum.clear();
    for (int n=0; n<nthrds; n++) {
      wgtSum.insert(wgtSum.end(), wgtVal[n]->begin(), wgtVal[n]->end());
    }

    // Only check for crazy collisions after first step
    //
    static bool notFirst = false;

    if (notFirst and numSanityStop) {

      MPI_Reduce(&Max, &gbMax,  1, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);

				// Too many collisions!!
      if (myid==0 && gbMax>numSanityMax) {
	std::cerr << printDivider << std::endl
		  << " *** Too many collisions in NTC: " << gbMax << std::endl
		  << printDivider << std::endl;
	raise(SIGTERM);		// Signal stop, please!
      }

    } else notFirst = true;

    MPI_Reduce(&Ovr, &accOvr, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Acc, &accAcc, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Tot, &accTot, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

    pHOT_iterator c(*tree);

    if (NTC_DIST) {

      while (c.nextCell()) {
	for (auto q : qv) {
	  try {
	    NTC::NTCitem::vcMap v = ntcdb[c.Cell()->mykey].CrsVel(q);
	    for (auto i : v) {
	      for (auto j : i.second)
		qq[i.first][j.first][q].push_back(j.second);
	    }
	  }
	  catch (NTC::NTCitem::Error &error) {}
	}
      }
      
      for (int n=1; n<numprocs; n++) {

	if (myid == n) {

	  unsigned sz = ntcSum.size();
	  MPI_Send(&sz, 1, MPI_UNSIGNED, 0, 224, MPI_COMM_WORLD);
	  if (sz) MPI_Send(&ntcSum[0], sz, MPI_DOUBLE, 0, 225, MPI_COMM_WORLD);
	  ntcSum.clear();

	  sz = wgtSum.size();
	  MPI_Send(&sz, 1, MPI_UNSIGNED, 0, 289, MPI_COMM_WORLD);
	  if (sz) MPI_Send(&wgtSum[0], sz, MPI_DOUBLE, 0, 290, MPI_COMM_WORLD);
	  wgtSum.clear();

	  sz = qq.size();
	  MPI_Send(&sz, 1, MPI_UNSIGNED, 0, 226, MPI_COMM_WORLD);
	  
	  for (auto i : qq) {

	    sKeyPair   k = i.first;
	    KeyConvert k1 (k.first );
	    KeyConvert k2 (k.second);
	    
	    int i1 = k1.getInt();
	    int i2 = k2.getInt();
	    int i3 = i.second.size();
	    
	    MPI_Send(&i1, 1, MPI_INT, 0, 227, MPI_COMM_WORLD);
	    MPI_Send(&i2, 1, MPI_INT, 0, 228, MPI_COMM_WORLD);
	    MPI_Send(&i3, 1, MPI_INT, 0, 229, MPI_COMM_WORLD);

	    int base = 230;
	    
	    for (auto j : i.second) {
	      NTC::T val = j.first;
	      unsigned num = j.second.size();
	      unsigned short uuu[5] = {std::get<0>(val),
				       std::get<1>(val).first,
				       std::get<1>(val).second,
				       std::get<2>(val).first,
				       std::get<2>(val).second};

	      MPI_Send(&uuu, 5, MPI_UNSIGNED_SHORT, 0, base+0, MPI_COMM_WORLD);
	      MPI_Send(&num, 1, MPI_UNSIGNED,       0, base+1, MPI_COMM_WORLD);

	      double z;		// temporary

	      for (auto s : j.second) {
		MPI_Send(&(z=s.first),    1, MPI_DOUBLE,   0, base+2, MPI_COMM_WORLD);
		unsigned num2 = s.second.size();
		MPI_Send(&num2,           1, MPI_UNSIGNED, 0, base+3, MPI_COMM_WORLD);
		if (num2) 
		  MPI_Send(&s.second[0], num2, MPI_DOUBLE,   0, base+4, MPI_COMM_WORLD);
	      }
	      base += 5;
	    }
	  }
	  
	} // END: process send to root
	
	if (myid==0) {
	  
	  std::vector<double> v;
	  unsigned short uuu[5];
	  unsigned sz, num, num2;
	  double dval;
	  int i1, i2, i3;
	
	  MPI_Recv(&sz, 1, MPI_UNSIGNED, n, 224, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  if (sz) {
	    v.resize(sz);
	    MPI_Recv(&v[0], sz, MPI_DOUBLE, n, 225, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    ntcSum.insert(ntcSum.end(), v.begin(), v.end());
	  }


	  MPI_Recv(&sz, 1, MPI_UNSIGNED, n, 289, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  if (sz) {
	    v.resize(sz);
	    MPI_Recv(&v[0], sz, MPI_DOUBLE, n, 290, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    wgtSum.insert(wgtSum.end(), v.begin(), v.end());
	  }

	  MPI_Recv(&sz, 1, MPI_UNSIGNED, n, 226, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  for (unsigned z=0; z<sz; z++) {

	    MPI_Recv(&i1, 1, MPI_INT, n, 227, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(&i2, 1, MPI_INT, n, 228, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(&i3, 1, MPI_INT, n, 229, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    
	    KeyConvert k1(i1);
	    KeyConvert k2(i2);
	    sKeyPair   k (k1.getKey(), k2.getKey());

	    int base = 230;

	    for (int w=0; w<i3; w++) {

	      MPI_Recv(&uuu, 5, MPI_UNSIGNED_SHORT, n, base+0, MPI_COMM_WORLD, 
		       MPI_STATUS_IGNORE);

	      MPI_Recv(&num, 1, MPI_UNSIGNED, n, base+1, MPI_COMM_WORLD, 
		       MPI_STATUS_IGNORE);

	      NTC::T utup {uuu[0], speciesKey(uuu[2], uuu[3]), speciesKey(uuu[5], uuu[6])};

	      if (num) {

		for (unsigned j=0; j<num; j++) {

		  MPI_Recv(&dval, 1, MPI_DOUBLE, n, base+2, MPI_COMM_WORLD, 
			   MPI_STATUS_IGNORE);

		  MPI_Recv(&num2, 1, MPI_UNSIGNED, n, base+3, MPI_COMM_WORLD, 
			   MPI_STATUS_IGNORE);

		  if (num2) {
		    v.resize(num2);
		    MPI_Recv(&v[0], num2, MPI_DOUBLE, n, base+4, MPI_COMM_WORLD,
			     MPI_STATUS_IGNORE);
		    
		    qq[k][utup][dval].insert( qq[k][utup][dval].end(), v.begin(), v.end() );
		  }
		}
	      }

	      base += 5;

	    } // Loop over quantiles

	  } // Loop over species

	} // Root receive loop

	MPI_Barrier(MPI_COMM_WORLD);

      } // Process loop

      if (myid==0) {

	if (wgtSum.size()) {
	  ntcHisto = ahistoDPtr(new AsciiHisto<double>(ntcSum, 20, 0.01, true));
	  ntcSum.clear();
	}

	if (wgtSum.size()) {
	  wgtHisto = ahistoDPtr(new AsciiHisto<double>(wgtSum, 20, 0.01, true));
	  wgtSum.clear();
	}

	for (auto &u : qq) {
	  for (auto &v : u.second) {
	    for (auto &z : v.second) {
	      qqHisto[u.first][v.first][z.first] = 
		ahistoDPtr(new AsciiHisto<double>(z.second, 20, 0.01));
	    }
	  }
	}
      }	
    }
  }
}

void Collide::NTCstanza(std::ostream& out, CrsVelMap& vals, int j,
			const std::string& lab, 
			const std::vector<double>& pcent)

{
  // Loop through the percentile list
  //
  for (auto p : pcent) {
    std::ostringstream sout;
    // Label the quantile
    //
    sout << lab << "(" << std::round(p*100.0) << ")";
    out << std::setw(18) << sout.str();
    // For each species pair, print the quantile
    //
    for (auto u : vals) {
      for (auto v : u.second) {
	size_t indx = static_cast<size_t>(std::floor(v.second[j].size()*p));
	out << std::scientific << std::setprecision(4)
	    << std::setw(12) << v.second[j][indx];
      }
    }
    out << std::endl;
  }
}

void Collide::NTCstats(std::ostream& out)
{
  const std::vector<double> pcent = {0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95};

  if (myid==0 and NTC and DEBUG_NTC) {

    size_t spc = 53;

    out << std::string(spc, '-') << std::endl
	<< "[NTC diagnostics]"   << std::endl << std::scientific << std::left
	<< std::setw(14) << " Time"     << std::setw(16) << tnow   << std::endl
	<< std::setw(14) << " Over"     << std::setw(16) << accOvr << std::endl
	<< std::setw(14) << " Accepted" << std::setw(16) << accAcc << std::endl
	<< std::setw(14) << " Total"    << std::setw(16) << accTot << std::endl
	<< std::fixed;

    if (accTot>0)
      out << std::setw(14) << " Ratio"    
	  << std::setw(16) << static_cast<double>(accAcc)/accTot << std::endl
	  << std::setw(14) << " Fail"     
	  << std::setw(16) << static_cast<double>(accOvr)/accTot << std::endl;

    out << std::string(spc, '-') << std::endl << std::right;

    if (ntcHisto) {
      out << "Full NTC distribution" << std::endl;
      (*ntcHisto)(out);
      out << std::string(spc, '-') << std::endl << std::right;
    }
      
    if (wgtHisto) {
      out << "Subspecies weight distribution" << std::endl;
      (*wgtHisto)(out);
      out << std::string(spc, '-') << std::endl << std::right;
    }

    if (qqHisto.size() > 0) {

      for (auto i : qqHisto) {
	sKeyPair   k  = i.first;
	speciesKey k1 = k.first;
	speciesKey k2 = k.second;
	
	std::ostringstream sout;
	if (k1.second==0 and k2.second==0)
	  sout << " Species: <" << k1.first << "|" << k2.first << ">";
	else
	  sout << " Species: <" << k1.first << "," << k1.second 
	       << "|" << k2.first << "," << k2.second << ">";

	for (auto v : i.second) {
	  
	  for (auto j : v.second) {

	    if (j.second.get()) {
	      out << std::endl << std::string(spc, '-') << std::endl
		  << std::left << std::fixed << sout.str() << std::endl
		  << " Interact: (" << labels[std::get<0>(v.first)]
		  << ", [" << std::setw(3) << std::get<1>(v.first).first << ", " << std::setw(3) << std::get<1>(v.first).second << "]"
		  << ", [" << std::setw(3) << std::get<1>(v.first).first << ", " << std::setw(3) << std::get<1>(v.first).second << "]"
		  << ")" << std::endl
		  << " Quantile: " << j.first << std::endl
		  << std::string(spc, '-') << std::endl
		  << std::left << std::scientific;
	      (*j.second)(out);
	      out << std::endl;
	    }
	  }
	}
      }
    }
  }
}


void Collide::mfpCLGather()
{
  // This vector will accumulate from all threads
  //
  std::vector<double> data;

  for (auto &v : mfpCLdata) {
    // Append info from all threads
    //
    data.insert(data.end(), v.begin(), v.end());
    // Remove all data for next step(s)
    //
    v.clear();
  }

  for (int i=1; i<numprocs; i++) {

    if (i == myid) {
      unsigned dNum = data.size();
      MPI_Send(&dNum,       1, MPI_UNSIGNED, 0, 435, MPI_COMM_WORLD);
      MPI_Send(&data[0], dNum, MPI_DOUBLE,   0, 436, MPI_COMM_WORLD);
    }
				// Root receives from Node i
    if (0 == myid) {
      unsigned dNum;
      MPI_Recv(&dNum,       1, MPI_UNSIGNED, i, 435, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      std::vector<double> vTmp(dNum);
      MPI_Recv(&vTmp[0], dNum, MPI_DOUBLE,   i, 436, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      data.insert(data.end(), vTmp.begin(), vTmp.end());
    }
    
    if (myid==0 and data.size()) {
      mfpclHist  = ahistoDPtr(new AsciiHisto<double>(data, 20, 0.01));
    }
  }
}
  

void Collide::mfpCLPrint(std::ostream& out)
{
  if (mfpclHist.get()) {

    // Print the header for mfpcl quantiles
    //
    out << std::endl << std::string(53, '-')  << std::endl
	<< "-----Cell length / MFP distribution------------------" 
	<< std::endl
	<< std::string(53, '-') << std::endl << std::left;
    (*mfpclHist)(out);
    out << std::string(53, '-')  << std::endl << std::endl;
  }
}


// Return random 3d unit vector
std::vector<double> Collide::unit_vector()
{
  enum UV {trig, gaus};		// Method choice
  static const UV uv(gaus);	// Pick the method
  std::vector<double> ret(3);	// Return vector

  if (uv == trig) {
    double cos_th = 1.0 - 2.0*unit(random_gen);
    double sin_th = sqrt(fabs(1.0 - cos_th*cos_th));
    double phi    = 2.0*M_PI*unit(random_gen);
    ret[0] = sin_th*cos(phi);
    ret[1] = sin_th*sin(phi);
    ret[2] = cos_th;
  } else {
    double nrm = 0.0;
    for (auto & v : ret) {
      v = norm(random_gen);
      nrm += v*v;
    }
    nrm = sqrt(nrm);
    for (auto & v : ret) v /= nrm;
  } 

  return ret;
}
