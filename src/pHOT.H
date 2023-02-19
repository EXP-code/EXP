#ifndef _pHOT_H
#define _pHOT_H

#include <limits>

#include <localmpi.H>
#include <Particle.H>
#include <Component.H>
#include <ParticleFerry.H>
#include <pHOT_types.H>
#include <pCell.H>
#include <Timer.H>

class CellDiag
{
public:
  unsigned ntre, bcel, btot, mind, maxd;
  double   navg, nvar, bavg, bvar;

  CellDiag(unsigned nt=0) : ntre(nt),  bcel(0),   btot(0),
			    mind(std::numeric_limits<int>::max()), maxd(0),
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

				// Tree element lists
				// ------------------
  key_indx keybods;		// Link body key to body index
  key_cell frontier;		// Link cell key to cell ptr
  key_key  bodycell;		// Link body key to cell key

				// Out of bounds list
  set<indx_type> oob;

  pCell *root;
  
  double volume;

  //! Send contents of this cell to another process
  void sendCell(key_type key, int to, unsigned num);

  //! Receive the contents of a cell from another process
  void recvCell(int from, unsigned num);

  //! Get the cell key for the beginning of the key-sorted particle list
  key_type getHeadKey();

  //! Get the cell key for the end of the key-sorted particle list
  key_type getTailKey();

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
  void rrMerge(vector<key_wght>& in, vector<key_wght>& out);
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
  static bool use_weight;

  //! Frontier iterator
  key_cell::iterator fit;
  bool reset;
  unsigned total_cells;

  ParticleFerryPtr pf;

  //! For debugging
  string debugf;

  key_type key_min, key_max;

  unsigned min_cell, max_cell, adjcnt, d_val, d_pval;
  unsigned cntr_total, cntr_new_key, cntr_mine, cntr_not_mine, cntr_ship;
  vector<unsigned> chist;
  change_list change;
  vector<key_type> kbeg, kfin, loclist;
  unsigned sumstep, sumzero;

  vector<double> offset;

  unsigned n_xchange, m_xchange, numkeys, numk;
  Timer timer_keymake, timer_xchange, timer_convert, timer_overlap;
  Timer timer_prepare, timer_cupdate, timer_scatter, timer_repartn;
  Timer timer_tadjust, timer_cellcul, timer_keycomp, timer_keybods;
  Timer timer_waiton1, timer_waiton2, timer_keynewc, timer_keyoldc;
  Timer timer_waiton0, timer_keysort, timer_keygenr, timer_diagdbg;

  vector<float>    keymk3, exchg3, cnvrt3, tovlp3, prepr3, updat3;
  vector<float>    scatr3, reprt3, tadjt3, celcl3, keycm3, keybd3;
  vector<float>    wait03, wait13, wait23, keync3, keyoc3, barri3;
  vector<float>    diagd3, keyst3, keygn3;
  vector<unsigned> numk3;

  map<pCell*, unsigned> clevlst; // Map cell ptrs to cell level
  vector< set<pCell*> > clevels; // List of cells at each level

  //@{
  //! For CellLevelList
  vector<unsigned> Pcnt, Plev;
  unsigned Nlev;
  //@}
  
  template <typename T> 
  void getQuant(vector<T>& in, vector<T>& out);

  void bomb(const string& membername, const string& msg);

  //static bool use_weight;

  static unsigned klen;

  static void qtile_initialize(); // Quantile set up

#if HAVE_LIBCUDA==1

  void keyProcessCuda(std::vector<key_wght> keys);

  void sortCombineOne(std::vector<key_type>& one,
		      std::vector<key_type>& two,
		      std::vector<key_type>& comb);

  void parallelMergeOne(std::vector<key_type>& initl,
			std::vector<key_type>& final);

#endif

public:

  //! Effort value hysteresis (default: 0.25)
  static double hystrs;

  //! 3-vector of sides of rectangular prism (default: [2, 2, 2])
  static std::vector<double> sides;

  //! Origin of the prism (default: [1, 1, 1])
  //! E.g. the default selction defines a cube with coordinates [-1, 1]^3
  static std::vector<double> offst;

  //! Quantiles for reporting distribution of time per process
  static std::vector<unsigned> qtile;

  //! Number of quantiles
  static unsigned ntile;

  //! Use key subsampling to improve run time in the key partioning step
  //! (default: true)

  static bool     sub_sample;

  //! Print out sample cell diagnostics (default: false)
  static bool     samp_debug;

  //! For communication
  static MPI_Datatype CellDiagType;

  //@{
  //! Species info
  //int species;
  //set<int> spec_list;
  sKeySet spec_list;
  //@}

  //! Constructor
  //pHOT(Component* C, int species, set<int> spec_list);
  pHOT(Component* C, sKeySet spec_list);

  //! Destructor
  ~pHOT();

  //! Get the particle key.  Out of bounds (OOB) particles have key "0u"
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
  Particle * Body(unsigned k) const
  {
    PartMap::iterator ie = cc->particles.end();
    PartMap::iterator it = cc->particles.find(k);
    if (it==ie) throw "pHOT::Body: particle not found";

    return it->second.get();
  }


  //! Tree geometry
  void setSides(double x, double y, double z) 
  {
    sides = {x, y, z};
    volume = x*y*z;
  }

  void setOffset(double x, double y, double z) 
  { 
    offst  = {x, y, z};
    offset = offst;
  }

  //! Number of cells in tree
  unsigned Number() { return frontier.size(); }
  unsigned TotalNumber() { return total_cells; }
  unsigned CellCount(double pctl);

  enum ChangeFlag {CREATE, REMOVE, DELETE, RECOMP};

  set<pCell*>& CLevels(int M) {
    if (clevels.size()==0) makeCellLevelList();
    return clevels[M];
  }

  void makeCellLevelList();
  void adjustCellLevelList(unsigned);
  void gatherCellLevelList();
  void printCellLevelList(ostream& out, const string& msg);
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
  void checkSampleCells(const std::string&);

  //! Check that clevlst and frontier are consistent
  void checkCellLevelList(const std::string&);

  //! Check that clevlst and clevels are consistent
  void checkLevelLists(const std::string&);

  //! Check that bodies refer to the same cells as the body lists
  bool checkCellClevelSanity(const std::string& msg, unsigned mlevel);

  //! Check that cells in the multilevel lists are also on the master list
  bool checkCellClevel(const std::string& msg, unsigned mlevel);

  //! Check that every particle has an entry in the key-body list
  bool checkKeybods(const std::string&);

  //! Check that every particle has an entry in the key-body list (at
  //! this level and below)
  bool checkPartKeybods(const std::string& msg, unsigned mlev);

  //! Check that every particle has an entry in the key-body list,
  //! body-cell list, and that the cell is on the frontier
  bool checkKeybodsFrontier(const std::string& msg);

  //! Check that every particle has an entry in the key-cell list
  bool checkBodycell(const std::string&);

  //! Check that every body in a cell is on the frontier
  bool checkCellFrontier(const std::string&);

  //! Check that every body in the tree points to a body in the particle list
  bool checkParticles(ostream &out, const std::string& msg, bool pc=true);

  //! Check that bodies are not duplicated on the particle list
  bool checkDupes1(const std::string& msg);

  //! Check that bodies are not duplicated in multiple cells
  bool checkDupes2();

  //! Check that every cell in the level list is on the frontier
  bool checkFrontier(ostream &out, const std::string& msg);

  //! Get the total kinetic energy and dispersion
  double totalKE(double& KEtot, double& KEdsp);

  //! Get the total gass mass
  void totalMass(unsigned& Counts, double& Mass);

  //! Sanity count check for body indices
  void checkIndices();

  //! Print the frontier and summary statistics
  void dumpFrontier(std::ostream& out);

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

  //! Check effort weighting in current partition
  void checkEffort(unsigned mlevel);
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
	      vector<float>    &tadjust, vector<float>    &celcull,
	      vector<float>    &keycomp, vector<float>    &keybods,
	      vector<float>    &keysort, vector<float>    &keygenr,
	      vector<float>    &waiton0, vector<float>    &waiton1,
	      vector<float>    &waiton2, vector<float>    &keynewc,
	      vector<float>    &keyoldc, vector<float>    &treebar,
	      vector<float>    &diagdbg, vector<unsigned> &numk    );
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

  double Volume() { return volume/( (key_type)1u << (3*fit->second->level) ); }

  void KE(double &tot, double &dsp) { fit->second->KE(tot, dsp); }

  double Mass() { return fit->second->Mass(); }

  Particle *Body(vector<unsigned long>::iterator k) { return fit->second->Body(k); }

};

#endif
