#ifndef _pCell_H
#define _pCell_H

#include <vector>
#include <list>
#include <set>
#include <map>
#include <deque>

#include "pHOT_types.H"

using namespace std;

class pHOT;

class pCell
{
private:

  pHOT* tree;
  deque<double> CRMlist;
  double CRMsum;
  unsigned CRMnum;

public:

  //! Target collision bucket size
  static unsigned bucket;

  //! Target macroscopic bucket size
  static unsigned Bucket;

  //! Number of bits per dimension
  static unsigned nbits;

  //! Number of entries in running CRM average (default: 8)
  static unsigned CRMcnt;

  //! Sampling method
  static bool NTC;

  //! Constructors
  pCell(pHOT* tr);
  pCell(pCell* mom, unsigned id);

  //! Destructor
  ~pCell();

  int owner;
  key_type mykey;
  key_type mask;
  unsigned level;
  unsigned count;
  unsigned maxplev;
  double   effort;

  /**
     @param 0 mass
     @param 1 vel2x
     @param 2 vel2y
     @param 3 vel2z
     @param 4 velx
     @param 5 vely
     @param 6 velz
     @param 7 posx
     @param 8 posy
     @param 9 posz
  */
  vector<double> state;

  pCell* parent;
  pCell* sample;
  map<unsigned, pCell*> children;
  bool isLeaf;
  
  key_indx keys;
  set<unsigned> bods;

  //! Add a body, return the current frontier node
  pCell* Add(const key_pair&, change_list* change=0);

  //! Remove a body
  bool Remove(const key_pair&, change_list* change=0);

  //! Remove all bodies
  void RemoveAll();

  /**
     Find the cell node in the current frontier that should contain 
     this body and return null if not found
  */
  pCell* findNode(const key_type&);

  //! Check that this key belongs to this branch
  bool isMine(const key_type& key) {
    if (key==0) return false;
    key_type sig = (key_type)(key - mask) >> 3*(nbits-level);
    if (sig!=0) return false;
    return true;
  }

  //! Update keys list
  void UpdateKeys(const key_pair& oldpair, const key_pair& newpair);

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

  //! Return the length scale factor from the quad tree
  double Scale();

  //! Mean position
  //@{
  void MeanPos(double &x, double &y, double& z);
  void MeanPos(vector<double>& p);
  //@}

  //! Mean velocity
  //@{
  void MeanVel(double &x, double &y, double& z);
  void MeanVel(vector<double>& v);
  //@}

  //! Remake the maxplev variable from the particles
  unsigned remake_plev();

  //! Find the cells whose body number is "just" larger than Bucket
  pCell* findSampleCell();

  //! Return number of bodies in Sample Cell
  unsigned SampleBodyNumber() { return sample->bods.size(); }

  //! Return bodies from this cell
  Particle* Body(set<unsigned>::iterator k);

  //! Maximum relative velocity average (<0 if no entries)
  double CRMavg();

  //! Add maximum to list
  void CRMadd(double crm);

};

#endif
