#ifndef _pHOT_H
#define _pHOT_H

#include <values.h>		// For MAXINT

#include <localmpi.h>
#include <Particle.H>
#include <Component.H>
#include <ParticleFerry.H>
#include <pHOT_types.H>
#include <pCell.H>

#include <Timer.h>
#include <BarrierWrapper.H>

class CellDiag
{
public:
  unsigned ntre, bcel, btot, mind, maxd;
  double   navg, nvar, bavg, bvar;

  CellDiag(unsigned nt=0) : ntre(nt),  bcel(0),   btot(0),
			    mind(MAXINT),         maxd(0),
			    navg(0.0), nvar(0.0), bavg(0.0), bvar(0.0)
  {}
};


class pHOT
{

  friend class pCell;
  friend class pHOT_iterator;

protected:

  Component *cc;

  unsigned number;

  key_indx keybods;
  key_cell frontier;
  key_key bodycell;

				// Out of bounds list
  set<indx_type> oob;

  pCell *root;
  
  double volume;

  void sendCell(key_type key, int to, unsigned num);
  void recvCell(int from, unsigned num);

  //@{
  /** Allow for adding multiple types of partitioning variants.  So
      far, I have implemented only the initial Hilbert curve
      partition.
  */
  enum PartitionType {Hilbert} partType;

  //! Standard Hibert curve partitioning algorithm
  void partitionKeysHilbert(vector<key_wght>& keys,
			    vector<key_type>& kbeg, vector<key_type>& kfin);

  //! Call the requested algorithm
  void partitionKeys(vector<key_wght>& keys,
		     vector<key_type>& kbeg, vector<key_type>& kfin) {
    if (partType==Hilbert)       partitionKeysHilbert (keys, kbeg, kfin);
    else                         partitionKeysHilbert (keys, kbeg, kfin);
  }
  //@}

  //@{
  //! For combining key-weight lists for partitioning
  void parallelMerge(vector<key_wght>& in, vector<key_wght>& out);
  void sortCombine(vector<key_wght>& one, vector<key_wght>& two, 
		   vector<key_wght>& comb);
  //@}

  //! Which processor owns the cell
  unsigned find_proc(vector<key_type>& keys, key_type key);

  //@{
  //! Partition the out-of-bounds particles among the procesors
  void spreadOOB();
  void checkOOB(vector<unsigned>&);
  //@}

  //@{
  //! State accumulation
  vector<int> numfront;
  vector<int> displace;
  //@}

  //! For repartition
  bool use_weight;

  //! Frontier iterator
  key_cell::iterator fit;
  bool reset;
  unsigned total_cells;

  ParticleFerry pf;

  //! For debugging
  key_type key_min, key_max;

  unsigned min_cell, max_cell, adjcnt, d_val, d_pval;
  unsigned cntr_total, cntr_new_key, cntr_mine, cntr_not_mine, cntr_ship;
  vector<unsigned> chist;
  change_list change;
  vector<key_type> kbeg, kfin, loclist;
  unsigned sumstep, sumzero;
  BarrierWrapper *barrier;

  double *offset;

  unsigned n_xchange, m_xchange, numkeys, numk;
  Timer timer_keymake, timer_xchange, timer_convert, timer_overlap;
  Timer timer_prepare, timer_cupdate, timer_scatter, timer_repartn;
  Timer timer_tadjust, timer_keycall, timer_keycomp, timer_keybods;
  Timer timer_waiton1, timer_waiton2, timer_keynewc, timer_keyoldc;
  Timer timer_waiton0, timer_diagdbg;

  vector<float>    keymk3, exchg3, cnvrt3, tovlp3, prepr3, updat3;
  vector<float>    scatr3, reprt3, tadjt3, keycl3, keycm3, keybd3;
  vector<float>    wait03, wait13, wait23, keync3, keyoc3, barri3;
  vector<float>    diagd3;
  vector<unsigned> numk3;

  template <typename T> 
  void getQuant(vector<T>& in, vector<T>& out);

  void bomb(const string& membername, const string& msg);

  static unsigned klen;

public:

  //! 3-vector of sides of rectangular prism (default: [2, 2, 2])
  static double   sides[];

  //! Origin of the prism (default: [1, 1, 1])
  //! E.g. the default selction defines a cube with coordinates [-1, 1]^3
  static double   offst[];

  //! Quantiles for reporting distribution of time per process
  static unsigned qtile[];

  //! Use key subsampling to improve run time in the key partioning step
  //! (default: true)

  static bool     sub_sample;

  //! Print out partition diagnostics (default: false)
  static bool     keys_debug;

  //! Print out sample cell diagnostics (default: false)
  static bool     samp_debug;

  //! For communication
  static MPI_Datatype CellDiagType;

  //! Constructor
  pHOT(Component* C);

  //! Destructor
  ~pHOT();

  key_type getKey(double *);
  static string printKey(key_type);

  //! Turn on/off weighted partitioning
  void setWeights(bool onoff) { use_weight = onoff; }

  /**
     Particle handling
  */
  //@{
  void makeTree();
  void Repartition(unsigned);
  void adjustTree(unsigned);
  Particle *Body(unsigned k) { return &cc->particles[k]; }


  //! Tree geometry
  void setSides(double x, double y, double z) 
  {
    sides[0] = x; sides[1] = y; sides[2] = z; 
  }

  void setOffset(double x, double y, double z) 
  { 
    offst[0] = x; offst[1] = y; offst[2] = z; 
    for (int k=0; k<3; k++) offset[k] = offst[k];
  }

  //! Number of cells in tree
  unsigned Number() { return frontier.size(); }
  unsigned TotalNumber() { return total_cells; }
  unsigned CellCount(double pctl);

  enum ChangeFlag {CREATE, REMOVE, KILL, RECOMP};
  map<pCell*, unsigned> clevlst;
  vector< set<pCell*> > clevels;
  void makeCellLevelList();
  void adjustCellLevelList(unsigned);
  void printCellLevelList(ostream& out);
  void computeCellStates();

  /**
     Density routines
  */
  //@{
  //! Do the extra work to determine the kinetic state per cell
  void makeState();

  //! Evaluate the kinetic state for the given input position
  void State(double *pos, double& dens, double& temp, 
	     double& vx, double& vy, double& vz);
  double minVol();
  double maxVol();
  double medianVol();
  double Volume() { return volume; }
  //@}


  /**
     Debugging and diagostics
  */
  //@{
  //! Print the distribution of cell levels on the frontier
  void logFrontierStats();
  //! Check that sample cells exist for all frontier cells
  void checkSampleCells(const char *);
  //! Check that clevlst and frontier are consistent
  void checkCellLevelList(const char *);
  //! Check that every particle has an entry in the key-body list
  bool checkKeybods();
  //! Check that every particle has an entry in the key-body list (at
  //! this level and below)
  bool checkPartKeybods(unsigned mlev);
  //! Check that every particle has an entry in the key-cell list
  bool checkBodycell();
  //! Check that every body in the tree points to a body in the particle list
  bool checkParticles(ostream &out, bool pc=true);
  //! Check that bodies are not duplicated in multiple cells
  bool checkDupes();
  //! Check that every cell in the level list is on the frontier
  bool checkFrontier(ostream &out);
  //! Get the total kinetic energy and dispersion
  double totalKE(double& KEtot, double& KEdsp);
  //! Get the total gass mass
  void totalMass(unsigned& Counts, double& Mass);
  //! Sanity count check for body indices
  void checkIndices();
  //! Print the frontier and summary statistics
  void dumpFrontier();
  //! 
  void densEmit(unsigned lev, pCell *p);
  void densCheck();
  void statFrontier();
  void countFrontier(vector<unsigned>& ncells, vector<unsigned>& bodies);
  unsigned checkNumber();
  //! Return the cached number of out-of-bounds particles
  unsigned oobNumber();
  void testFrontier(string& filename);
  void Slice(int nx, int ny, int nz, string cut, string prefix);
  void Slab(vector<int>& n, vector<double>& pmin, vector<double>& pmax,
	    string cut, 
	    vector<double>&    x, vector<double>& dens, vector<double>& temp, 
	    vector<double>& velx, vector<double>& vely, vector<double>& velz);
  //! Scan for and count out-of-bounds particles
  void checkBounds(double, const char *);
  void adjustCounts(ostream& out);
  double checkAdjust() {
    double result = 0.0;
    if (myid) return result;
    if (sumstep) {
      result = static_cast<double>(sumzero)/sumstep;
      sumstep = sumzero = 0;
    }
    return result;
  }
  void Xchange(unsigned int& nX, unsigned int& mX) { 
    nX = n_xchange;
    mX = m_xchange;
    n_xchange = m_xchange = 0;
  }

  bool onFrontier(key_type k) {
    if (frontier.find(k) == frontier.end()) return false;
    else return true;
  }

  //@}


  //! Return timing info for essential routines
  //@{
  //! Collect timing from all processes
  void CollectTiming();

  //! Process and return timing stats from all processes
  void Timing(vector<float>    &keymake, vector<float>    &exchange, 
	      vector<float>    &convert, vector<float>    &overlap, 
	      vector<float>    &prepare, vector<float>    &update,
	      vector<float>    &scatter, vector<float>    &repartn,
	      vector<float>    &tadjust, vector<float>    &keycall,
	      vector<float>    &keycomp, vector<float>    &keybods,
	      vector<float>    &waiton0, vector<float>    &waiton1,
	      vector<float>    &waiton2, vector<float>    &keynewc,
	      vector<float>    &keyoldc, vector<float>    &treebar, 
	      vector<float>    &diagdbg, vector<unsigned> &numk);
  //@}

};

/**
   Reentrant iterator
*/
class pHOT_iterator
{
private:
  //! Calling tree
  pHOT *tr;

  //! Frontier iterator
  key_cell::iterator fit;

  //! Volume of tree
  double volume;
  
  bool first;

public:

  pHOT_iterator(pHOT& p) : 
    tr(&p), fit(p.frontier.begin()), volume(p.volume), first(true) {}

  unsigned nextCell() 
  { 
    if (first) first = false;
    else       fit++;
    if (fit==tr->frontier.end()) return 0;
    return fit->second->bods.size(); 
  }

  pCell* Cell() { 
    if (fit==tr->frontier.end()) {
#ifdef DEBUG
      cout << "No cell!" << endl;
#endif
      return 0;
    }
    return fit->second; 
  }

#ifdef INT128
  double Volume() { return volume/(( (key_type)1u << (3*fit->second->level) ).toDouble()); }
#else
  double Volume() { return volume/( (key_type)1u << (3*fit->second->level) ); }
#endif
  void KE(double &tot, double &dsp) { fit->second->KE(tot, dsp); }

  double Mass() { return fit->second->Mass(); }

  Particle *Body(vector<unsigned long>::iterator k) { return fit->second->Body(k); }

};

#endif
