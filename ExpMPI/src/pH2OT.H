#ifndef _pH2OT_H
#define _pH2OT_H

#include <localmpi.h>
#include <Particle.H>
#include <Component.H>
#include <ParticleFerry.H>
#include <pHOT_types.H>
#include <pCell.H>

#include <Timer.h>
#include <BarrierWrapper.H>

#include <values.h>

typedef pair<unsigned, int> pkey_proc;

class pH2OT;
class pTree;
class tTree;

class tCell : public sCell
{
  friend class tTree;
  friend class pTree;

protected:

  static unsigned pbits;
  static tTree *tr;
  void getRange(vector<double>& min3, vector<double>& max3);

public:

  //! Track number of instances
  static int live;

  //! Construtors
  tCell(tTree* tree, unsigned nb);
  tCell(tCell* mom,  unsigned id);

  //! Destructor
  ~tCell();

  //! Walk up to find the cell with at least @param bucket entries
  sCell *findSampleCell(unsigned bucket);

  //! Find the cell node in the current frontier that should contain
  //! this point and return null if not found
  tCell* findNode(const double* p);

  //! Find the cell node in the current frontier that should contain
  //! this key and return null if not found
  tCell* findNode(const unsigned& key);

  //! Find the cell node in the current frontier that should contain
  //! this point and add it if necessary
  tCell* findAddNode(const double* p);

  //! Find the cell node in the current frontier that should contain
  //! this key and add it if necessary
  tCell* findAddNode(const unsigned& key);

  //! Add state vector to tree, accumulated from leave to root
  void Accum(const double* s);

  //! Get the cell key containing point *p
  static unsigned getKey(const double *p);

  vector<double> getCorner();

  //! Public data
  //@{
  unsigned mykey;
  unsigned mask;
  pTree*   ptree;
  tCell*   parent;
  map<unsigned, tCell*> children;
  //@}

  //! Return total volume of cell
  double Volume();

  //! Return the length scale factor from the quad tree
  double Scale();

  //! Check that this key belongs to this branch
  bool isMine(const unsigned& key) {
    if (key==0u) return false;
    unsigned sig = (key - mask) >> 3*(pbits-level);
    if (sig!=0u) return false;
    return true;
  }

  //! To which child does this belong
  unsigned childId(unsigned key)
  {
    unsigned id = key - mask;
    id >>= 3*(pbits - 1 - level);
    if (id>7) {
      cerr << "Crazy value: " << id << endl;
    }
    return id;
  }
};

class tTree
{
  friend class tCell;
  friend class pH2OT;

protected:

  pH2OT *ph;
  unsigned maxlevel, oob;
  tCell *root;
  void printTree(tCell *c);
  tCell* TreeWalk(tCell *c, unsigned key);
  
public:
  //! Track number of instances
  static int live;

  tTree() : maxlevel(0), oob(0), root(0) { live++; }
  tTree(pH2OT* p, unsigned pbits);
  ~tTree() { live--; delete root; }

  map<unsigned, tCell*> frontier;

  void reset(pH2OT* p, unsigned pbits);
  void ZeroAccum() { ZeroAccum(root); };
  void ZeroAccum(tCell *t);
  void AddState(const double* s);
  void AddCell(tCell *cell);
  void printTree() { printTree(root); }
  void printCell(const double* p);
  unsigned getOOB() { unsigned ret=oob; oob=0; return ret; }

  tCell* MakeCell(unsigned key);

  pTree* RetrievePtree(unsigned key) {
    if (frontier.find(key) == frontier.end()) {
      frontier[key] = MakeCell(key);
    }
    return frontier[key]->ptree;
  }

};


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

class pTree
{
  friend class tCell;
  friend class pCell;
  friend class pH2OT;

protected:

  pH2OT* ph;
  tCell* pc;
  unsigned tkey;
  map<pCell*, unsigned> clevlst;
  vector<unsigned> chist;
  change_list change;

public:
  //! Track number of instances
  static int live;

  pCell *root;
  vector<double> C;
  key_indx keybods;
  key_cell frontier;
  key_key  bodycell;
  unsigned minc, maxc;
  vector< set<pCell*> > clevels;
  enum ChangeFlag {CREATE, REMOVE, KILL, RECOMP};

  //! Constructor
  pTree(pH2OT* p, tCell *q, unsigned key);

  //! Destructor
  ~pTree();

  //! Main members
  void makeTree();
  void computeStatus(unsigned mlevel, vector< vector<unsigned long> >& exchange);
  void cellUpdate(unsigned mlevel);
  void makeCellLevelList();
  void adjustCellLevelList(unsigned mlevel);
  unsigned getCHist(unsigned n) {
    if (n<minc || n>maxc) return 0u;
    return chist[n - minc];
  }

  //! Diagnostics
  void densEmit(unsigned lev, pCell *p);
  bool checkBodycell();
  bool checkLevelList();
  void getRange(vector<double>& x, vector<double>& X) { pc->getRange(x, X); }
};


class pH2OT
{
  friend class pTree;
  friend class pCell;
  friend class pH2OT_iterator;

protected:

  Component *cc;
  unsigned number;
  set<unsigned long> oob;

  unsigned int pkmask, tkmask, tottrees;
  double pkfactor;
  set<pkey_proc> lproc;

  map<unsigned, int> keyproc;
  void Remap_by_number();
  void Remap_by_effort();
  template <typename T> 
  void Remap_diag(map<unsigned, T>& count);

  int getProc(unsigned key);
  vector<double> getCorner(unsigned);
  
  double total_volume, volume;

  //! For repartition
  bool use_weight, use_number;
  void particleExchange(vector< vector<unsigned long> >& exchange);

  //! Frontier iterator
  key_cell::iterator fit;
  unsigned total_cells;

  //! State accumulation
  //@{
  vector<unsigned> numfront;
  //@}

  ParticleFerry pf;

  //! For debugging
  key_type key_min, key_max;

  unsigned min_cell, max_cell, d_val, d_pval;
  vector<unsigned> chist;
  BarrierWrapper *barrier;

  unsigned numkeys;
  Timer timer_cstatus, timer_xchange, timer_convert, timer_tadjust;
  Timer timer_prepare, timer_cupdate, timer_scatter, timer_repartn;
  Timer timer_keybods, timer_schecks, timer_waiton0, timer_waiton1;
  Timer timer_waiton2, timer_bodlist, timer_celladj, timer_getsta1;
  Timer timer_getsta2, timer_getsta3;
  
  vector<float>    cstat3, exchg3, prepr3, cnvrt3, tadjt3, updat3;
  vector<float>    scatr3, reprt3, keybd3, schks3, wait03, wait13;
  vector<float>    wait23, bodls3, cladj3, barri3, cells1, cells2;
  vector<float>    cells3;
  vector<unsigned> numk3;

  template <typename T> 
  void getQuant(vector<T>& in, vector<T>& out);

  //! Stuff to be done when the global volume is changed
  void geomUpdate();
  //! Remake and accumulate the cell statistics
  void computeCellStates();

  void bomb(const string& membername, const string& msg);

  static unsigned klen;

public:

  //! 3-vector of sides of rectangular prism (default: [2, 2, 2])
  static double   sides[];

  //! 3-vector of sides of oct-tree (default: [2, 2, 2])
  static double   box[];

  //! Origin of the prism (default: [1, 1, 1])
  //! E.g. the default selction defines a cube with coordinates [-1, 1]^3
  static double   offset[];

  //! Use key subsampling to improve run time in the key partioning step
  //! (default: true)
  static bool     sub_sample;

  //! Quantiles for reporting distribution of time per process
  static unsigned qtile[];

  //! Do consistency check diagnostics (default: false)
  static bool     list_check;

  //! For communication
  static MPI_Datatype CellDiagType;


  //! Constructor
  pH2OT(Component* C);

  //! Destructor
  ~pH2OT();

  key_type getKey(double *, vector<double>&);
  static string printKey(key_type);

  //! Turn on/off weighted partitioning
  void setWeights(bool onoff) { use_weight = onoff; }

  /**
     Particle handling
  */
  //@{
  tTree trees;

  void Remap();
  void makeTree();
  void adjustTree(unsigned);
  void Repartition(unsigned);

  Particle *Body(unsigned k) { return &cc->particles[k]; }


  //! Tree geometry
  inline void setSides(double x, double y, double z) 
  {
    sides[0] = x; sides[1] = y; sides[2] = z; 
				// Reset total volume
    total_volume = sides[0]*sides[1]*sides[2];
				// Reset subbox size
    for (int k=0; k<3; k++) box[k] = sides[k]*pkfactor;
				// Reset subbox volume
    volume = box[0]*box[1]*box[2];
				// Reset subbox corners
    geomUpdate();
  }

  inline void setOffset(double x, double y, double z) 
  { 
    offset[0] = x; offset[1] = y; offset[2] = z; 
				// Reset subbox corners
    geomUpdate();
  }

  //! Number of cells in tree
  inline unsigned Number() { 
    unsigned n=0;
    for (map<unsigned, tCell*>::iterator 
	   t=trees.frontier.begin();  t!=trees.frontier.end(); t++) 
      n += t->second->ptree->frontier.size();
    return n;
  }

  //! Return the total number of cells
  unsigned TotalNumber() { return total_cells; }

  //! Return the percentile @param pctl number of bodies in the cells
  unsigned CellCount(double pctl);

  //! Distribute out-of-bounds bodies among all processors
  void spreadOOB();

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

  //! Return the volume of the top level HOT cell; that is the root
  //! cell in each subtree
  double Volume() { return volume; }

  //! Return the volume of the entire structure
  double TotalVolume() { return total_volume; }
  //@}


  /**
     Debugging and diagostics
  */
  //@{
  //! Return the total kinetic energy for all particles in all cells
  double totalKE(double& KEtot, double& KEdsp);

  //! Return the total gas mass for all particles in all cells
  void totalMass(unsigned& Count, double& Mass);

  //! Compute the density for cells
  void densCheck();

  //! Return the number cells and bodies on the frontier
  void countFrontier(vector<unsigned>& ncells, vector<unsigned>& bodies);

  //! Return the number of particles in all cells
  unsigned checkNumber();

  //! Return the number of out-of-bounds particles
  unsigned oobNumber();

  //! Return summary statistics for frontier cells in @param filename
  void testFrontier(string& filename);

  //! Check that all of the cells assigned to levels are on the
  //! frontier
  bool checkFrontier(ostream &out);

  
  //! Compute summary statistics for all cells
  void checkCellTree();

  //! Check that all of the particles are accounted for in cells
  //@{
  bool checkParticles(ostream &out, bool pc=true);
  bool checkParticles() { return checkParticles(cout); }
  //@}
  
  void checkBounds(double, const char *);

  //! Report on out-of-bounds particles
  void checkOOB(vector<unsigned>& sendlist);

  //! Check that all the bodies in cells in the level lists are in the
  //! bodycell list
  bool checkLevelList();

  //! Check that all the particles with non-zero keys belong to a
  //! bodycell list
  bool checkBodycell();

  //! Print the cell level list counts
  void printCellLevelList(ostream &out);

  //! Return timing info for essential routines
  //@{
  //! Collect timing from all processes
  void CollectTiming();

  //! Process and return timing stats from all processes
  void Timing(vector<float> &cstatus, vector<float> &exchange, 
	      vector<float> &prepare, vector<float> &convert, 
	      vector<float> &tadjust, vector<float> &update,
	      vector<float> &scatter, vector<float> &repartn,  
	      vector<float> &keybods, vector<float> &schecks, 
	      vector<float> &waiton0, vector<float> &waiton1,
	      vector<float> &waiton2, vector<float> &bodlist,  
	      vector<float> &celladj, vector<float> &getsta1, 
	      vector<float> &getsta2, vector<float> &getsta3,
	      vector<float> &barrier, vector<unsigned> &numk);
  //@}

  //! Toggle the diagostic barrier checking (off for production)
  inline void MPIchk(bool tf) { if (tf) barrier->on(); else barrier->off(); }

  //@}
};

/**
   Reentrant iterator
*/
class pH2OT_iterator
{
private:
  //! Calling tree
  pH2OT *pr;

  //! Tree iterator
  map<unsigned, tCell*>::iterator tr;

  //! Frontier iterator
  key_cell::iterator fit;

  //! First iteration
  bool first;

  //! Volume of tree
  double volume;

public:

  pH2OT_iterator(pH2OT& p) : pr(&p), tr(p.trees.frontier.begin()), first(true), volume(p.volume) {}

  unsigned nextCell();

  pCell* Cell() { return fit->second; }
  
#ifdef I128
  double Volume() { return volume/(( (key_type)1u << (3*fit->second->level) ).toDouble()); }
#else
  double Volume() { return volume/( (key_type)1u << (3*fit->second->level) ); }
#endif

  void KE(double &tot, double &dsp) { fit->second->KE(tot, dsp); }

  double Mass() { return fit->second->Mass(); }

  Particle *Body(set<unsigned long>::iterator k) { return fit->second->Body(k); }

};

#endif
