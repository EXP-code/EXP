#ifndef _pCell_H
#define _pCell_H

#include <vector>
#include <map>

#include "pHOT_types.H"

using namespace std;

class pHOT;

class pCell
{
private:

  pHOT* tree;

public:

  //! Target collision bucket size
  static unsigned bucket;

  //! Target macroscopic bucket size
  static unsigned Bucket;

  //! Number of bits per dimension
  static unsigned nbits;

  //! Constructors
  pCell(pHOT* tr);
  pCell(pCell* mom, unsigned id);

  //! Destructor
  ~pCell();

  int owner;
  key_type mykey;
  key_type mask;
  double volume;
  unsigned level;
  unsigned count;

  /**
     @param 0 mass
     @param 1 vel2
     @param 2 velx
     @param 3 vely
     @param 4 velz
  */
  vector<double> state;

  pCell* parent;
  pCell* sample;
  map<unsigned, pCell*> children;
  bool isLeaf;
  
  key_set keys;
  vector<unsigned> bods;

  //! Add a body, return the current frontier node
  pCell* Add(const key_pair&);

  //! Recursively zero the tree's state vector
  void zeroState();
  
  //! Walk down, accumulating state
  void accumState();
  void accumState(unsigned count, vector<double>& _state);

  //! Compute state values for the maximum level
  void Find(key_type key, unsigned& curcnt, unsigned& lev,
	    vector<double>& state);

  //! Compute the node key at the child level
  unsigned childId(key_type key);

  //! Return kinetic energy per unit mass per cell (total & temperature)
  void KE(double &tot, double &dsp);

  //! Return velocity statistics
  void Vel(double &mass, vector<double>& v, vector<double>& v2);

  //! Return total mass in cell;
  double Mass();

  //! Return total volume of cell
  double Volume();

  //! Find the cells whose body number is "just" larger than Bucket
  pCell* findSampleCell();

  //! Return number of bodies in Sample Cell
  unsigned SampleBodyNumber() { return sample->bods.size(); }

  //! Return bodies from Sample Cell
  Particle* SampleBody(unsigned int k);

  //! Return bodies from this cell
  Particle* Body(unsigned int k);

};

#endif
