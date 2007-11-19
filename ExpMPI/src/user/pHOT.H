#ifndef _pHOT_H
#define _pHOT_H

#include "localmpi.h"
#include "Particle.H"
#include "ParticleFerry.H"
#include "pHOT_types.H"
#include "pCell.H"

class pHOT
{

  friend class pCell;
  friend class pHOT_iterator;

protected:

  unsigned number;

  key_indx keybods;
  key_cell frontier;
  key_key bodycell;
  indx_part bodies;

  vector<Particle> *data;
  pCell *root;
  
  double volume;

  void sendCell(key_type key, int to, unsigned num);
  void recvCell(int from, unsigned num);

  //! Frontier iterator
  key_cell::iterator fit;
  bool reset;

  ParticleFerry pf;

  //! For debugging
  key_type key_min, key_max;

public:

  static unsigned neg_half;
  static double s[];

  pHOT();
  ~pHOT();

  key_type getKey(double *);
  static string printKey(key_type);

  /**
     Particle handling
  */
  //@{
  void makeTree();
  void Repartition();
  void sendParticles(vector<Particle>& p);
  void gatherParticles(vector<Particle>& p);
  Partstruct *Body(unsigned k) { return &bodies[k]; }

  //! Number of cells in tree
  unsigned Number() { return frontier.size(); }

  /**
     Density routines
  */
  //@{
  //! Do the extra work to determine the kinetic state per cell
  void makeState();

  //! Evaluate the kinetic state for the given input position
  void State(double *pos, double& dens, double& temp, 
	     double& vx, double& vx, double& vz);
  double minVol();
  double medianVol();
  double Volume() { return volume; }
  //@}


  /**
     Debugging and diagostics
  */
  //@{
  bool checkKeybods();
  void dumpFrontier();
  void densEmit(unsigned lev, pCell *p);
  void densCheck();
  void statFrontier();
  void testFrontier(string& filename);
  void testDump();
  void Slice(int nx, int ny, int nz, string cut, string prefix);
  void Slab(vector<int>& n, vector<double>& pmin, vector<double>& pmax,
	    string cut, 
	    vector<double>&    x, vector<double>& dens, vector<double>& temp, 
	    vector<double>& velx, vector<double>& vely, vector<double>& velz);
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
    return fit->second->bods.size(); 
  }

  pCell* Cell() { return fit->second; }

  double Volume() { return volume/(1 << (3*fit->second->level)); }
  void KE(double &tot, double &dsp) { fit->second->KE(tot, dsp); }

  double Mass() { return fit->second->Mass(); }

  Partstruct *Body(unsigned k) { return &(tr->bodies[fit->second->bods[k]]); }
};

#endif
