#ifndef _Collide_H
#define _Collide_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include <boost/tuple/tuple.hpp>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include "Particle.H"
#include "pHOT.H"
#include "Timer.h"
#include "NTC.H"

class Collide;


//! Helper class to pass info to threaded member
struct thrd_pass_arguments
{
  //! Pointer to iterator
  pHOT *tree;

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
  //! Pointer to iterator
  pHOT *tree;

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
   The collision base class for DSMC
   
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
*/
class Collide
{
  
public:
  
  //! Diagnostic record type
  typedef std::pair<double, vector<double> > Precord;

  //! Collide return tuple
  typedef boost::tuple<unsigned, unsigned>   UU;
  

protected:

  //! Fiducial component
  Component *c0;

  //! Calling force
  ExternalForce *caller;
  
  //@{
  //! Threading stuff
  int nthrds;
  thrd_pass_Collide *td;
  thrd_pass_tstep   *tdT;
  pthread_t *t;
  //@}
  
  //@{
  //! Pseudorandom generation
  ACG *gen;
  Uniform *unit;
  Normal *norm;
  //@}
  
  //! Cell diagnostic records indexed by thread
  vector<Precord> prec;
  
  //! collide return counters
  UU ret;

  //@{
  // Timing and debug
  typedef std::list< std::pair<long, unsigned> > effortL;

  unsigned            ncells, bodycount, stepcount;
  Timer               forkTime, snglTime, waitTime, joinTime, diagTime;
  TimeElapsed         forkSoFar, snglSoFar, waitSoFar;
  TimeElapsed         diagSoFar, joinSoFar;
  vector<double>      forkSum, snglSum, waitSum, diagSum, joinSum;
  bool                effortAccum;
  vector<effortL>     effortNumber;
  vector<unsigned>    ntcOvr, ntcTot;
  unsigned            accOvr, accTot;

  boost::tuple<double, double, double> ntcRes;
  
  //@}

  //! Quantile index function
  unsigned Qi(double q, unsigned sz) {
    if (q>=1.0) return sz-1;
    if (q<=0.0) return 0;
    return (unsigned)floor(q*sz);
  }
  
  //@{
  //! For timing the fraction spent in threaded force methods

  vector<struct timeval> tv_list;
  vector<double> timer_list;
  
  inline void thread_timing_beg(int id)
  {
    if (VERBOSE>5) {
      gettimeofday(&tv_list[id], 0);
      timer_list[2*id] = 
	tv_list[id].tv_usec*1.0e-6 +
	(tv_list[id].tv_sec % 1000);
    }
  }
  
  inline void thread_timing_end(int id)
  {
    if (VERBOSE>5) {
      gettimeofday(&tv_list[id], 0);
      timer_list[2*id+1] = 
	tv_list[id].tv_usec*1.0e-6 +
	(tv_list[id].tv_sec % 1000);
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
  
  std::vector< vector<pCell*> > cellist;
  
  std::vector<double> disptot;
  double masstot;
  double excessCtot, excessEtot;

  //@}

  
  //@{
  //! Parameters set by caller and internally

  double   diamfac;
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
  
  //! Initializes cross sections
  virtual void initialize_cell
  (pHOT* const tree, pCell* const cell, double vrelmax, int id) = 0;
  
  //! Additional DSMC-specific initialization
  virtual void initialize_cell_dsmc
  (pHOT* const tree, pCell* const cell, 
   sKey2Umap& nsel, double vrelmax, double tau, int id) = 0;
  
  //! Additional EPSM-specific initialization
  virtual void initialize_cell_epsm
  (pHOT* const tree, pCell* const cell, sKey2Umap& nsel, 
   double vrelmax, double tau, int id) = 0;
  
  //! Compute energy adjustments diagnostics for the cell after
  //! collisions/cooling
  virtual void finalize_cell(pHOT* const tree, pCell* const cell, 
			     double kedsp, int id) {}
  
  //! Return total scattering + elastic cross section "matrix" for
  //! species-species interactions
  virtual sKey2Dmap& totalCrossSections(int id) = 0;

  //! Return the full cross section for particle selection computation
  virtual double crossSection(pHOT* const tree, 
			      Particle* const p1, Particle* const p2,
			      double cr, int id=0) { return 0; }
  
  //! Return the particle selection matrix
  virtual sKey2Umap generateSelection
  (pCell* const c, sKeyDmap* const Fn, double crm, double tau, int id,
   double& meanLambda, double& meanCollP, double& totalNsel) = 0;

  //! Do the inelastic computation
  virtual int inelastic(pHOT* const tree, 
			Particle* const p1, Particle* const p2,
			double *cr, int id=0) { return 0; }
  
  virtual void 
  velocityUpdate(Particle* const p1, Particle* const p2, double cr);

  //! Has LTE heating and cooling algorithm
  virtual bool hasHeatCool() = 0;
  
  //! Compute the EPSM algorithm for the given cell
  virtual void EPSM(pHOT* const tree, pCell* const cell, int id);
  
  //! Use EPSM if true
  bool use_epsm;
  
  //@{
  //! For diagnostic output
  vector<double> lostSoFar_EPSM;
  
  void debug_list(pHOT& tree);
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
  vector<TimeElapsed> initSoFar, collSoFar, elasSoFar, cellSoFar, epsmSoFar;
  vector<TimeElapsed> stat1SoFar, stat2SoFar, stat3SoFar, coolSoFar, listSoFar;
  vector<TimeElapsed> curcSoFar;
  vector<int>         collCnt, numbSum;
  vector<double>      initSum, collSum, elasSum, cellSum, epsmSum;
  vector<double>      stat1Sum, stat2Sum, stat3Sum, coolSum, listSum;
  vector<long>        minUsage, maxUsage, minPart, maxPart, minCollP, maxCollP;
  
  static const int nEPSMT = 6;
  vector< vector<Timer> > EPSMT;
  vector< vector<TimeElapsed> > EPSMTSoFar;
  vector<long> EPSMtime, CPUH;
  //@}
  
  //! Species map output
  std::string species_file_debug;

  //! Cached molecular weight (computed once; could be done per cell
  //! but temperature is a mean quantity anyway)
  double mol_weight;

  //! Initialize atomic weights (public data)
  void atomic_weights_init();


  //! NTC database
  NTCdb ntcdb;

public:
  
  //! Integrate the cooling solution
  static bool ESOL;		
  
  //! Use the original Pullin velocity selection algorithm
  static bool PULLIN;		
  
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
  
  /** Use EPSM when mean free path/cell size ratio is below EPSMratio
      (default = -1.0, never) */
  static double EPSMratio;
  
  //! Minimum number of particles for EPSM
  static unsigned EPSMmin;
  
  //! Temperature floor in resampling dispersion for EPSM
  static double TFLOOR;
  
  //! Artificially increase or decrease cooling rate (default: 1, no change)
  static double ENHANCE;
  
  //! CrossSection-Velocity factor for NTC database
  static double NTCAVG;

  //! Atomic weights (public data)
  static std::map<unsigned short, double> atomic_weights;

  /** Constructor
      @param force is the caller
      @param comp is the fiducial component
      @param hsdiam is the hard-sphere diameter in Bohr units
      @param diamfac is the diameter scaling
  */
  Collide(ExternalForce *force, Component *comp, 
	  double hsdiam, double diamfac, int nth=1);
  
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
  //! Collision thread
  void *collide_thread(void *arg);
  void collide_thread_fork(pHOT *tree, sKeyDmap* Fn);
  //@}
  
  //@{
  //! Timestep thread
  virtual void *timestep_thread(void *arg);
  void compute_timestep(pHOT* tree, double coolfrac);
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
			    double T=0);
  //@}
  
  //@{ 
  //! Collision diagnostics. Currently, these are only in CollideIon
  virtual void printCollGather()  {}
  virtual void printCollSummary() {}
  virtual void resetColls()       {}
  //@}

  //@{ 
  //! Auxilliary diagnostics.
  virtual void auxGather() {}
  virtual void auxPrint(ostream& out) {}
  //@}

  //@{
  //! Cell distribution of average (cross*vel)
  std::map< sKeyPair, std::vector<double> > both, cros, crat;
  //}

  //! Gather stats for No-Time-Counter algorithm
  virtual void NTCgather(pHOT* const tree);

  //! Print stats for No-Time-Counter algorithm
  virtual void NTCstats(ostream& out);

protected:

  //! Helper function for NTCstats output
  void NTCstanza(std::ostream& out, 
		 std::map< sKeyPair, std::vector<double> >& vals,
		 const std::string& lab, const std::vector<double>& pcent);
};

#endif
