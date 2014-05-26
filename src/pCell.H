#ifndef _pCell_H
#define _pCell_H

#include <vector>
#include <list>
#include <set>
#include <map>
#include <deque>

#include <global.H>
#include <pHOT_types.H>

using namespace std;

class pHOT;
class Component;

template< class T1, class T2 >
std::pair<T1, T2> operator+(const pair<T1, T2>& a, const pair<T1, T2>& b)
{
  std::pair<T1, T2> ret(a.first + b.first, a.second + b.second);
  return ret;
};

template< class T1, class T2 >
std::pair<T1, T2> operator-(const pair<T1, T2>& a, const pair<T1, T2>& b)
{
  std::pair<T1, T2> ret(a.first - b.first, a.second - b.second);
  return ret;
};

template< class T1, class T2 >
std::pair<T1, T2> operator+=(pair<T1, T2>& a, const pair<T1, T2>& b)
{
  a.first  += b.first;
  a.second += b.second;
  return a;
};

template< class T1, class T2 >
std::pair<T1, T2> operator-=(pair<T1, T2>& a, const pair<T1, T2>& b)
{
  a.first  -= b.first;
  a.second -= b.second;
  return a;
};

class sCell
{
public:
  typedef std::pair<double, double> dPair;

protected:
  typedef std::deque<dPair>         ddqPair;

  ddqPair       VelCrsList;
  dPair         VelCrsSum;
  unsigned      VelCrsNum;
  
  std::map<speciesKey, ddqPair>     VelCrsListM;
  std::map<speciesKey, dPair>       VelCrsSumM;
  std::map<speciesKey, unsigned>    VelCrsNumM;
  
public:
  
  // Maximum number of cached relative velocities
  static size_t VelCrsSZ; 
  
  int owner;
  unsigned level;
  std::map<speciesKey, unsigned> count;
  unsigned ctotal;
  
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
  //map<int, vector<double> > state;
  map<speciesKey,vector<double> > state;
  vector<double> stotal;
  
  //@{
  //! Maximum relative velocity average (<0 if no entries)
  dPair VelCrsAvg(speciesKey indx);
  dPair VelCrsAvg();
  //@}
  
  //@{
  //! Add maximum to list
  void VelCrsAdd(speciesKey indx, const dPair& val);
  void VelCrsAdd(const dPair& val);
  void VelCrsAdd(speciesKey indx, double vrm, double crm) 
  { VelCrsAdd(indx, dPair(vrm, crm)); }
  void VelCrsAdd(double vrm, double crm)
  { VelCrsAdd(dPair(vrm, crm)); }
  //@}
  
  //@{
  //! Return kinetic energy per unit mass per cell (total & temperature)
  void KE(speciesKey indx, double &tot, double &dsp);
  void KE(double &tot, double &dsp);
  //@}
  
  //@{
  //! Return total mass in cell
  double Mass();
  double Mass(speciesKey);
  //@}
  
  //@{
  //! Return total number in cell
  unsigned Count();
  unsigned Count(speciesKey);
  //@}
  
  //! Return total volume of cell
  virtual double Volume() = 0;
  
  //! Return the length scale factor from the quad tree
  virtual double Scale() = 0;
  
  //! Mean position
  //@{
  void MeanPos(speciesKey indx, double &x, double &y, double& z);
  void MeanPos(speciesKey indx, vector<double>& p);
  void MeanPos(double &x, double &y, double& z);
  void MeanPos(vector<double>& p);
  //@}
  
  //! Mean velocity
  //@{
  void MeanVel(speciesKey indx, double &x, double &y, double& z);
  void MeanVel(speciesKey indx, vector<double>& v);
  void MeanVel(double &x, double &y, double& z);
  void MeanVel(vector<double>& v);
  //@}
};

class pCell : public sCell
{
private:
  
  pHOT* tree;
  Component *C;
  
public:
  
  //! Track number of instances
  static int live;
  
  //! Target collision bucket size
  static unsigned bucket;
  
  //! Target macroscopic bucket size
  static unsigned Bucket;
  
  //! Maximum # of cell expansions to find macroscopic sample cell
  static unsigned deltaL;
  
  //! Sampling method
  static bool NTC;
  
  //! Constructors
  pCell(pHOT* tr);
  pCell(pCell* mom, unsigned id);
  
  //! Destructor
  ~pCell();
  
  key_type mykey;
  key_type mask;
  unsigned maxplev;
  
  pCell* parent;
  sCell* sample;
  map<unsigned, pCell*> children;
  bool isLeaf;
  
  key_indx keys;
  std::vector<unsigned long> bods;
  
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
    if (key==0u) return false;
    key_type sig = (key_type)(key - mask) >> 3*(nbits-level);
    if (sig!=0u) return false;
    return true;
  }
  
  //! Remove key from list
  void RemoveKey(const key_pair& oldpair);
  
  //! Update keys list
  void UpdateKeys(const key_pair& oldpair, const key_pair& newpair);
  
  //! Recursively zero the tree's state vector
  void zeroState();
  
  //! Walk down, accumulating state
  void accumState();
  void accumState(sKeyUmap& count, sKeyvDmap& _state);
  
  //@{
  //! Compute state values for the maximum level
  void Find(key_type key, unsigned& curcnt, unsigned& lev,
	    vector<double>& state);
  
  void Find(key_type key, sKeyUmap& curcnt, unsigned& lev,
	    sKeyvDmap& st);
  //@}
  
  
  //! Compute the node key at the child level
  unsigned childId(key_type key);
  
  //! Remake the maxplev variable from the particles
  unsigned remake_plev();
  
  //! Find the cells whose body number is "just" larger than Bucket
  sCell* findSampleCell();
  
  //! Return bodies from this cell
  Particle* Body(vector<unsigned long>::iterator k);
  
  //! Return velocity statistics for particles in the cell
  void Vel(double &mass, vector<double>& v, vector<double>& v2);
  
  //! Return total volume of cell
  double Volume();
  
  //! Return the length scale factor from the quad tree
  double Scale();
  
};

const speciesKey defaultKey(0, 0);

#endif
