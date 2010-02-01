#ifndef _pHOT_H
#define _pHOT_H

#include <localmpi.h>
#include <Particle.H>
#include <Component.H>
#include <ParticleFerry.H>
#include <pHOT_types.H>
#include <pCell.H>

#include <Timer.h>
#include <BarrierWrapper.H>

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

  void partitionKeys(vector<key_wght>& keys,
		     vector<key_type>& kbeg, vector<key_type>& kfin);
  void parallelMerge(vector<key_wght>& in, vector<key_wght>& out);
  void sortCombine(vector<key_wght>& one, vector<key_wght>& two, 
		   vector<key_wght>& comb);

  unsigned find_proc(vector<key_type>& keys, key_type key);
  void spreadOOB();
  void checkOOB(vector<unsigned>&);

  //! For repartition
  bool use_weight;

  //! Frontier iterator
  key_cell::iterator fit;
  bool reset;
  unsigned total_cells;

  ParticleFerry pf;

  //! For debugging
  key_type key_min, key_max;

  unsigned min_cell, max_cell, adjcnt;
  vector<unsigned> chist;
  change_list change;
  vector<key_type> kbeg, kfin, loclist;
  unsigned sumstep, sumzero;
  BarrierWrapper *barrier;

  ACG *gen;
  Uniform *unit;
  double *offset;


  unsigned n_xchange, m_xchange, numkeys;
  Timer timer_keymake, timer_xchange, timer_convert, timer_overlap;
  Timer timer_prepare, timer_cupdate, timer_scatter, timer_repartn;
  Timer timer_tadjust, timer_keycall, timer_keycomp, timer_keybods;
  Timer timer_waiton1, timer_waiton2;

  void bomb(const string& membername, const string& msg);

  static unsigned klen;

public:

  //! 3-vector of sides of rectangular prism (default: [2, 2, 2])
  static double   sides[];
  //! Origin of the prism (default: [1, 1, 1])
  //! E.g. the default selction defines a cube with coordinates [-1, 1]^3
  static double   offst[];
  //! Additional randomly selected offset with given extent to prevent
  //! gridding effects (default: none)
  static double   jittr[];
  //! Use key subsampling to improve run time in the key partioning step
  //! (default: true)
  static bool     sub_sample;
  //! Print out partition diagnostics (default: false)
  static bool     keys_debug;

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


  //! Number of cells in tree
  unsigned Number() { return frontier.size(); }
  unsigned TotalNumber() { return total_cells; }
  unsigned CellCount(double pctl);

  enum ChangeFlag {CREATE, REMOVE, KILL, RECOMP};
  map<pCell*, unsigned> clevlst;
  vector< set<pCell*> > clevels;
  void makeCellLevelList();
  void adjustCellLevelList(unsigned);

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
  bool checkKeybods();
  bool checkPartKeybods(unsigned mlev);
  bool checkBodycell();
  bool checkParticles(bool pc=true);
  bool checkDupes();
  bool checkFrontier();
  double totalKE(double& KEtot, double& KEdsp);
  void totalMass(unsigned& Counts, double& Mass);
  void checkIndices();
  void dumpFrontier();
  void densEmit(unsigned lev, pCell *p);
  void densCheck();
  void statFrontier();
  void countFrontier(vector<unsigned>& ncells, vector<unsigned>& bodies);
  unsigned checkNumber();
  unsigned oobNumber();
  void testFrontier(string& filename);
  void Slice(int nx, int ny, int nz, string cut, string prefix);
  void Slab(vector<int>& n, vector<double>& pmin, vector<double>& pmax,
	    string cut, 
	    vector<double>&    x, vector<double>& dens, vector<double>& temp, 
	    vector<double>& velx, vector<double>& vely, vector<double>& velz);
  void checkBounds(double, const char *);
  void adjustTiming(double &keymake, double &xchange, double &prepare, 
		    double &convert, double &overlap, double &update,
		    double &scatter, double &repartn, double &tadjust,
		    double &keycall, double &keycomp, double &keybods,
		    double &waiton1, double &waiton2, unsigned &numkeys);
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
    if (fit==tr->frontier.end()) {
#ifdef DEBUG
      cout << "End of cell list!" << endl;
#endif
      return 0;
    }
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

  double Volume() { return volume/((key_type)1 << (3*fit->second->level)); }
  void KE(double &tot, double &dsp) { fit->second->KE(tot, dsp); }

  double Mass() { return fit->second->Mass(); }

  Particle *Body(set<unsigned>::iterator k) { return fit->second->Body(k); }

};

#endif
