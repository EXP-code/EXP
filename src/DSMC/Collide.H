#ifndef _Collide_H
#define _Collide_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <tuple>
#include <map>


#include "ExternalForce.H"
#include "AsciiHisto.H"
#include "Particle.H"
#include "pHOT.H"
#include "Timer.H"
#include "NTC.H"

#include <circular_buffer.H>

class Collide;


//! Helper class to pass info to threaded member
struct thrd_pass_arguments
{
  //! Effective number of particles
  sKeyDmap *fn;

  //! Thread counter id
  int id;
};

struct thrd_pass_Collide 
{
  //! Pointer to parent
  Collide *p;

  //! Arguments
  thrd_pass_arguments arg;
};


//! Helper class to pass info to threaded member
struct tstep_pass_arguments
{
  //! Cooling fraction
  double coolfrac;

  //! Thread counter id
  int id;
};

struct thrd_pass_tstep
{
  //! Pointer to parent
  Collide *p;

  //! Arguments
  tstep_pass_arguments arg;
};


/**
   \brief The collision base class for DSMC
   
   <h2> Virtual functions that must be implemented by a derived class </h2>
   <ul>
   
   <li> initialize_cell_dsmc().  This function initializes cross
   sections (and possibly) heating and cooling rates for use in DSMC.
   
   <li> initialize_cell_epsm(). Initializes mean heating and cooling
   rates for use in EPSM.
   
   <li> finalize_cell(). Computes diagnostics for the cell after
   collisions/cooling.
   
   <li> totalCrossSections(). Return total cross section "matrix" for
   species-species interactions.  The "matrix" is implemented as an
   STL map<int, map<int, double> > so sparse indices do not matter but
   we can still use the matrix-like syntatic construction
   matrix[i1][i2].
   
   <li> crossSection(). Return the cross section for the collision for
   use in the particle selection algorithm.

   <li> inelastic().  Performs the inelastic computation.  Returns 0 on
   success.
   
   <li> hasHeatCool().  Return true if the derived class can provide a
   heating and cooling function for EPSM.  Otherwise, return false.
   
   </ul>

   <h2> Requirements for DSMC </h2>
   <ul>

   <li> The no-time-counter (NTC), Monte-Carlo-based collision
   selector assumes that the <cross section>x<velocity> value is the
   maximum possible in the cell.  This value is maintained per cell by
   the NTCdb class.  Evaluations are entered into the database
   sequentially and trimmed on accorded to FIFO.

   <li> The Collide class provides diagnostics for the NTC algorithm
   that are printed in the <runtag>.DSMC_log file.

   <li> Time steps to be smaller than the mean time of flight through
   the cell.

   <li> The collision cells are current chosen to have a target number
   of particles per cell (based on current DSMC best practice, this is
   chosen to be of order 10).  

   <li> The number of particles should be chosen so that the ratio of
   the mean-free-path (MFP) for collisions to collision cell length
   scale is of order unity (typically MFP/Length approx. 0.3 is
   suggested).  This requires choosing an appropriate number of
   particles.

   <li> However, other trade-offs may be more appropriate and should
   be implemented.  It is possible that MFP/Length maybe be smaller or
   larger than 0.3 for the target number of particles per cell.  So,
   one may wish to choose cell size so that MPF/Length < 0.3 even if
   the particle count becomes larger than the target.

   </ul>
   
*/
class Collide
{
  
protected:

  std::string name_id;
  std::string version_id;

public:
  
  //! Diagnostic record type
  typedef std::pair<double, vector<double> > Precord;

  //! Collide return tuple
  typedef std::tuple<unsigned, unsigned> UU;
  
  //@{
  //! Electron species key
  const static unsigned short electronTyp;
  const speciesKey electronDef = {electronTyp, electronTyp};
  const sKeyPair   electronKey = {electronDef, electronDef};
  //@}

  //@{
  //! Maximum number of collisions per step

  //! Allow counter to stop job
  static bool numSanityStop;

  //! Maximum number per step
  static unsigned numSanityMax;

  //! Verbose messaging
  static bool numSanityMsg;

  //! Lower thresh for reporting
  static unsigned numSanityVal;

  //! Upper thresh for reporting
  static unsigned numSanityFreq;
  
  //! Random seed for initialization
  static unsigned seed;

  //! Toggle for deep debugging output for hybrid method (set to false
  //! for production)
  static constexpr bool init_dbg = false;

  //! Initial time for deep debugging output
  static constexpr double init_dbg_time = 0.0;

  //! Species selection for deep debugging output
  static constexpr unsigned short init_dbg_Z = 2;

  //! Using cuda (default: false)
  static bool use_cuda;

  //@}

  //! Interaction-type labels
  class Labels : public std::map<int, std::string>
  {
  public:
    static std::string def;

    std::string& operator[] (int i)
    {
      std::map<int, std::string> *p = 
	static_cast<std::map<int, std::string>*>(this);
      if (p->find(i) == p->end()) (*p)[i] = def;
      return (*p)[i];
    }

  } labels;

protected:

  //! Fiducial component
  Component *c0;

  //! Calling force
  ExternalForce *caller;
  
  //! Tree pointer cached from component
  pHOT *tree;

  //! Called with multistep level
  int mlev;

  //@{
  //! Threading stuff
  int nthrds;
  thrd_pass_Collide *td;
  thrd_pass_tstep   *tdT;
  pthread_t *t;
  pthread_mutex_t tlock;
  //@}
  
  //@{
  //! Pseudorandom generation
  std::uniform_real_distribution<> unit;
  std::normal_distribution<>       norm;
  //@}
  
  //! Cell diagnostic records indexed by thread
  std::vector<Precord> prec;
  
  //! collide return counters
  UU ret;

  //@{
  // Timing and debug
  typedef std::list< std::pair<double, unsigned> > effortL;

  unsigned            ncells, bodycount, stepcount;
  Timer               forkTime, snglTime, waitTime, joinTime, diagTime;
  double              forkSoFar, snglSoFar, waitSoFar;
  double              diagSoFar, joinSoFar;
  vector<double>      forkSum, snglSum, waitSum, diagSum, joinSum;
  bool                effortAccum;
  vector<effortL>     effortNumber;
  vector<unsigned>    ntcOvr, ntcAcc, ntcTot;
  unsigned            accOvr, accAcc, accTot;

  std::tuple<double, double, double> ntcRes;
  
  //@}

  //! Quantile index function
  unsigned Qi(double q, unsigned sz) {
    if (q>=1.0) return sz-1;
    if (q<=0.0) return 0;
    return (unsigned)floor(q*sz);
  }
  
  //@{
  //! For timing the fraction spent in threaded force methods

  vector<std::time_t> timer_list;
  
  inline void thread_timing_beg(int id)
  {
    if (VERBOSE>5) {
      timer_list[2*id] = std::chrono::system_clock::to_time_t(std::chrono::high_resolution_clock::now());
    }
  }
  
  inline void thread_timing_end(int id)
  {
    if (VERBOSE>5) {
      timer_list[2*id+1] = std::chrono::system_clock::to_time_t(std::chrono::high_resolution_clock::now());
    }
  }
  //@}
  
  //@{
  //! Diagnostic variables

  std::vector<unsigned> numcnt;
  std::vector<unsigned> selcnt;
  std::vector<unsigned> colcnt;
  std::vector<double>   tsrat, kerat, derat, tdens, tvolc, ttemp, tdelt, tseln;
  
  static const unsigned Nphase = 4, Nmfp = 5;
  std::vector<Precord> tphase, tmfpst;
  std::vector< vector<Precord> > tphaseT, tmfpstT;
  
  std::vector< vector<unsigned> > numcntT;
  std::vector< vector<unsigned> > colcntT;
  std::vector<unsigned> error1T;
  std::vector<unsigned> sel1T, col1T, epsm1T, Nepsm1T;
  std::vector<double> KEtotT, KElostT, tmassT, decolT, decelT, exesCT, exesET;
  std::vector< vector<double> > mfpratT, tsratT, tdensT, tvolcT;
  std::vector< vector<double> > ttempT, tselnT, keratT, deratT;
  std::vector< vector<double> > tdeltT, tdispT, EoverT, CoverT;
  const static unsigned nvold = 6, numdiag = 11;
  std::vector<unsigned> tdiag, tdiag1, tdiag0;
  std::vector<unsigned> tcool, tcool1, tcool0;
  std::vector<double>   Eover, Eover1, Eover0;
  std::vector<double>   Cover, Cover1, Cover0;
  std::vector< vector<unsigned> > tdiagT, tcoolT;
  std::vector<unsigned> Vcnt, Vcnt1, Vcnt0;
  std::vector< vector<unsigned> > VcntT;
  std::vector<double>   Vdbl, Vdbl1, Vdbl0;
  std::vector< vector<double> > VdblT;

  std::vector<NTC::Interact<double>> retCrs;

  std::vector< vector<pCell*> > cellist;
  
  std::map<pCell*, double> selMFP;

  std::vector<double> disptot;
  double masstot;
  double excessCtot, excessEtot;
  //@}

  //@{
  //! For histogram of mfpCL values
  static bool MFPCL;
  std::vector< std::vector<double> > mfpCLdata;
  ahistoDPtr mfpclHist;
  //@}
  
  //@{
  //! Parameters set by caller and internally

  double   crossfac;
  double   hsdiam;
  unsigned seltot;
  unsigned coltot;
  unsigned epsmtot;
  unsigned epsmcells;
  unsigned errtot;
  unsigned numtot;
  //@}
  
  //@{
  //! Runtime features set by caller
  int use_key, use_temp, use_dens, use_delt, use_Kn, use_St, use_exes;
  int use_Eint;
  //@}
  
  //@{
  //! Diagnostic for entire cell iteration loop
  virtual void pre_cell_loop(int id)  {}
  virtual void post_cell_loop(int id) {}
  //@}

  //! Initializes cross sections
  virtual void initialize_cell
  (pCell* const cell, double vrelmax, double tau, int id) = 0;
  
  //! Additional DSMC-specific initialization
  virtual void initialize_cell_dsmc
  (pCell* const cell, 
   NTC::InteractD& nsel, double vrelmax, double tau, int id) = 0;
  
  //! Additional EPSM-specific initialization
  virtual void initialize_cell_epsm
  (pCell* const cell, NTC::InteractD& nsel, 
   double vrelmax, double tau, int id) = 0;
  
  //! Compute energy adjustments diagnostics for the cell after
  //! collisions/cooling
  virtual void finalize_cell(pCell* const cell, 
			     sKeyDmap* const Fn, double kedsp, double tau, 
			     int id) {}
  
  //! Return total scattering + elastic cross section "matrix" for
  //! species-species interactions
  virtual NTC::InteractD& totalCrossSections(pCell* const c, double cr, int id) = 0;

  //! Return the full cross section for particle selection computation
  virtual double crossSection(int id, pCell* const c, 
			      Particle* const p1, Particle* const p2,
			      double cr, const NTC::T& ityp=NTC::single) = 0;
  
  //! Initialize quantities for cross-section computation
  virtual void
  pairInfo(int id, pCell* const c,
	   Particle* const p1, Particle* const p2, double cr) = 0;


  //! Select an interaction pair
  virtual std::pair<int, int> pairSelect(int id, NTC::T T, int nbods)
  {
    std::pair<int, int> ret;
    ret.first  = std::min<int>(int(floor(unit(random_gen)* nbods   )), nbods-1);
    ret.second = std::min<int>(int(floor(unit(random_gen)*(nbods-1))), nbods-2);
    if (ret.second >= ret.first) ret.second++;
    return ret;
  }

  //! Return the particle selection matrix
  virtual NTC::InteractD generateSelection
  (pCell* const c, double crm, sKeyDmap* const Fn, double tau,
   double& meanLambda, double& meanCollP, double& totalNsel, int id) = 0;

  //! Do the inelastic computation
  virtual int inelastic(int id, pCell* const cell, 
			Particle* const p1, Particle* const p2, double *cr, 
			double wght, const NTC::T& ityp=NTC::single) = 0;
  
  virtual void 
  velocityUpdate(Particle* const p1, Particle* const p2, double cr);

  //! Has LTE heating and cooling algorithm
  virtual bool hasHeatCool() = 0;
  
  //! Compute the EPSM algorithm for the given cell
  virtual void EPSM(pCell* const cell, int id);
  
  //! Use EPSM if true
  bool use_epsm;
  
  //@{
  //! For diagnostic output
  vector<double> lostSoFar_EPSM;
  
  void debug_list();
  virtual void list_sizes();
  void list_sizes_proc(ostream*);
  void pre_collide_diag();
  unsigned post_collide_diag();
  
  template<typename T>
  void colldeHelper(ostream& out, const char* lab, vector<T>& v);
  //@}
  
  //@{
  //! Timing and debug (access for derived classes)
  vector<Timer>       initTime, collTime, elasTime, cellTime, curcTime, epsmTime;
  vector<Timer>       stat1Time, stat2Time, stat3Time, coolTime, listTime;
  vector<double>      initSoFar, collSoFar, elasSoFar, cellSoFar, epsmSoFar;
  vector<double>      stat1SoFar, stat2SoFar, stat3SoFar, coolSoFar, listSoFar;
  vector<double>      curcSoFar;
  vector<int>         collCnt, numbSum;
  vector<double>      initSum, collSum, elasSum, cellSum, epsmSum;
  vector<double>      stat1Sum, stat2Sum, stat3Sum, coolSum, listSum;
  vector<long>        minUsage, maxUsage, minPart, maxPart, minCollP, maxCollP;
  
  static const int nEPSMT = 6;
  vector< vector<Timer>  > EPSMT;
  vector< vector<double> > EPSMTSoFar;
  vector<double> EPSMtime, CPUH;
  //@}
  
  //! Species map output
  std::string species_file_debug;

  //! Cached molecular weight (computed once; could be done per cell
  //! but temperature is a mean quantity anyway)
  double mol_weight;

  //! Initialize atomic weights (public data)
  void atomic_weights_init();

  //@{
  //! NTC database
  NTC::NTCdb ntcdb;
  bool use_ntcdb;
  //@}

  //@{
  //! NTC threshold
  static double ntcThreshDef;
  double ntcThresh, ntcFactor;
  //@}

  //@{
  //! NTC diagnostoics
  const size_t bufCap = 16384;
  typedef circular_buffer<double> circBuf;
  typedef std::shared_ptr<circBuf> circBufPtr;
  std::vector<circBufPtr> ntcVal, wgtVal;
  std::vector<double> wgtSum, ntcSum;
  //@}

  //! Return random 3d unit vector
  std::vector<double> unit_vector();

  //! Print divider
  static std::string printDivider;

public:
  
  //! Integrate the cooling solution
  static bool ESOL;		
  
  //! Use the original Pullin velocity selection algorithm
  static bool PULLIN;		

  //! Use sample cell for value of mean collision velocity
  static bool SampleCRM;

  //! Print out sorted cell parameters
  static bool SORTED;
  
  //! Print out extra collision statistics for debugging
  static bool EXTRA;
  
  //! Print out T-rho plane for cells will mass weighting
  static bool PHASE;		
  
  //! Turn off collisions (for testing, debugging, etc.)
  static bool DRYRUN;           
  
  //! Set cooling to zero for testing (default: false)
  static bool NOCOOL;
  
  /** Use ensemble excess method: share uncooled internal energy among
      all particles in a cell.  If false, compute this on a
      particle-by-particle basis. */
  static bool ENSEXES;
  
  //! Sample based on maximum (true) or estimate from variance (false);
  static bool NTC;	        
  
  //! Sample based on maximum (true) or estimate from variance (false)
  //! without database;
  static bool NTCnodb;	        
  
  //! Accumulate time step diagnostics (default: false)
  static bool TSDIAG;	        
  
  //! Accumulate cell-volume diagnostics (default: false)
  static bool VOLDIAG;
  
  //! Power of two interval for ke/cool ratio histogram (default: 4)
  static int TSPOW;	        
  
  //! Accumulate mean free path diagnostics (default: false)
  static bool MFPDIAG;
  
  //! Use a mean-free-path estimate of the local time step (default: false)
  static bool MFPTS;
  
  //! Augment effort value by work done in collisions (default: true)  
  static bool EFFORT;
  
  //! Print out additional CPU timing detail
  static bool TIMING;

  //! Enable verbose NTC diagnostics (default: false)
  static bool DEBUG_NTC;

  //! Enable full NTC distribution output (default: false)
  static bool NTC_DIST;

  /** Use EPSM when mean free path/cell size ratio is below EPSMratio
      (default = -1.0, never) */
  static double EPSMratio;
  
  //! Minimum number of particles for EPSM
  static unsigned EPSMmin;
  
  //! Temperature floor in resampling dispersion for EPSM
  static double TFLOOR;
  
  //! Artificially increase or decrease cooling rate (default: 1, no change)
  static double ENHANCE;
  
  //! Atomic weights (public data)
  static std::map<unsigned, double> atomic_weights;

  //! Target number of collisions per cell (for time step selection)
  static unsigned collTnum;

  //! Testspread over sampling factor (default: 40)
  static int TestSpreadCount;	        
  
  /** Constructor
      @param force is the caller
      @param comp is the fiducial component
      @param hsdiam is the hard-sphere diameter in Bohr units
      @param crossfac is the cross-section scaling for elastic scattering
      @param name_id is the class name for identification
      @param version_id is the version identification info
      @param nth is the number of desired threads
  */
  Collide(ExternalForce *force, Component *comp, 
	  double hsdiam, double crossfac,
	  const std::string& name_id,
	  const std::string& version_id,
	  int nth=1);
  
  //! Destructor
  virtual ~Collide();
  
  //! Hard-sphere diameter
  virtual double hsDiameter();

  virtual void
  set_key(int pos) { use_key=pos; }
  
  virtual void 
  set_temp_dens(int Tpos, int Dpos) { use_temp=Tpos; use_dens=Dpos; }
  
  virtual void
  set_excess(int pos) { use_exes=pos; }
  
  virtual void set_timestep(int pos) { }
  
  virtual void
  set_Kn(int pos) { use_Kn=pos; }
  
  virtual void
  set_St(int pos) { use_St=pos; }
  
  virtual void
  set_Eint(int pos) { use_Eint=pos; }
  
  virtual void
  set_MFPTS(bool flag) { MFPTS = flag; }
  
  virtual double getCoolingRate(int id) = 0;
  
  unsigned select() { unsigned t=seltot; seltot=0; return t; }
  unsigned total()  { unsigned t=coltot; coltot=0; return t; }
  unsigned EPSMtotal()  { unsigned t=epsmtot; epsmtot=0; return t; }
  unsigned EPSMcells()  { unsigned t=epsmcells; epsmcells=0; return t; }
  unsigned errors() { unsigned t=errtot; errtot=0; return t; }
  unsigned cellN()  { return numtot; }
  unsigned medianNumber();
  unsigned medianColl();

  //! Return the number of bodies at this level and the number of collisions
  const UU& collide(pHOT& cl, sKeyDmap& Fn, int mlev=0, bool diag=false);
  
  void collQuantile(vector<double>& quantiles, vector<double>& coll_);
  void mfpsizeQuantile(vector<double>& quantiles, vector<double>& mfp_, 
		       vector<double>& ts_,       vector<double>& coll_,
		       vector<double>& cool_,     vector<double>& rate_,
		       unsigned& collnum,         unsigned& coolnum     );
  
  void energyExcess(double& ExesColl, double& ExesEPSM);
  void dispersion(vector<double>& disp);
  
  virtual double Etotal() { cerr << "Oops [1]" << endl; return 0.0; }
  virtual double Mtotal() { cerr << "Oops [2]" << endl; return 0.0; }
  virtual void Elost(double* collide, double* epsm)  
  { cerr << "Oops [3]" << endl; *collide=0; *epsm=0; }

  //@{
  //! Multihreading of collisions

  //! Collision thread.  Virtual to accomodate CUDA implementation
  virtual void *collide_thread(void *arg);

  //! Launch the threads
  void collide_thread_fork(sKeyDmap* Fn);
  //@}
  
  //@{
  //! Timestep thread
  virtual void *timestep_thread(void *arg);
  void compute_timestep(double coolfrac);
  //@}
  
  //@{
  //! Debug timing routines
  vector<double> dgnoscTime() { return diagSum; }
  vector<double> threadTime() { return forkSum; }
  vector<double> joinedTime() { return snglSum; }
  vector<double> waitngTime() { return waitSum; }
  vector<double> joinngTime() { return joinSum; }

  void EPSMtimingGather();
  void CPUHogGather();
  void EPSMtiming(ostream& out);
  void CPUHog(ostream& out);
  virtual void KElossGather() {}
  virtual void KEloss(ostream& out) {}
  //@}
  
  //@{
  //! Debug output routines
  void CollectTiming();
  void colldeTime(ostream& out);
  void voldiag(ostream& out);
  void tsdiag(ostream& out);
  //@}

  /** Molecular weight in atomic mass units.  This could be computed
      per cell but now, it is computed once to start, possibly using the
      particle distribution.  
  */
  virtual double molWeight();
  
  //@{
  //! Print out species counts
  virtual void gatherSpecies() {}
  virtual void printSpecies(std::map<speciesKey, unsigned long>& spec,
			    double T = 0.0);
  //@}
  
  //@{ 
  //! Collision diagnostics. Currently, these are only in CollideIon
  virtual void printCollGather()  {}
  virtual void printCollSummary() {}
  virtual void resetColls()       {}
  //@}

  //@{
  //! Expensive diagnostics, only use these if you must
  void mfpCLGather();
  void mfpCLPrint(std::ostream&out );
  //@}

  //@{ 
  //! Auxilliary diagnostics.
  virtual void auxGather() {}
  virtual void auxPrint(std::ostream& out) {}
  //@}

  //@{
  //! Cell distribution of average (cross*vel)
  typedef std::map<NTC::T, std::map<float, ahistoDPtr> > NTCuuHist;
  typedef std::map<double, std::vector<double> >         dvMap;
  typedef std::map<sKeyPair, NTCuuHist>                  NTCqqHist;
  typedef std::map<NTC::T, dvMap>                        qqMap;
  typedef std::map<sKeyPair, qqMap>                      CrsVelMap;

  //! Desired quantiles
  const std::vector<double> qv = {0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99};

  //! Store ntcdb data
  CrsVelMap qq;

  //! Histogram of ntcdb data
  NTCqqHist  qqHisto;

  //! Histogram of all target ratios
  ahistoDPtr ntcHisto;

  //! Histogram for subspecies wghts
  ahistoDPtr wgtHisto;
  //}

  //! Gather stats for No-Time-Counter algorithm
  virtual void NTCgather();

  //! Print stats for No-Time-Counter algorithm
  virtual void NTCstats(ostream& out);

  //! Cache NTC database for restart
  virtual void finish() { ntcdb.finish(); }

  //! Age out NTC database
  virtual void ageout() { ntcdb.ageOut(); }

  //! Return particle map from the component
  PartMap& Particles() { return c0->Particles(); }

protected:

  //! Helper function for NTCstats output
  void NTCstanza(std::ostream& out, CrsVelMap& vals, int n,
		 const std::string& lab, const std::vector<double>& pcent);
};

#endif
