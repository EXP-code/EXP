#ifndef _pCell_H
#define _pCell_H

#include <vector>
#include <tuple>
#include <list>
#include <set>
#include <map>
#include <deque>

#include <global.H>
#include <pHOT_types.H>

using namespace std;

class pHOT;
class Component;

//
// Tuple arithmetic
//

//----------------------------------------------------------------------
// Add two tuples
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator+(const std::tuple<T1, T2, T3>& a, 
				 const std::tuple<T1, T2, T3>& b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) + std::get<0>(b), 
			      std::get<1>(a) + std::get<1>(b), 
			      std::get<2>(a) + std::get<2>(b) );
  return ret;
};

//----------------------------------------------------------------------
// Subtract two tuples
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator-(const std::tuple<T1, T2, T3>& a, 
				 const std::tuple<T1, T2, T3>& b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) - std::get<0>(b), 
			      std::get<1>(a) - std::get<1>(b), 
			      std::get<2>(a) - std::get<2>(b) );
  return ret;
};

//----------------------------------------------------------------------
// Multiply a tuple by a scalar
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S>
std::tuple<T1, T2, T3> operator*(const std::tuple<T1, T2, T3>& a, const S b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) * b,
			   std::get<1>(a) * b,
			   std::get<2>(a) * b );
  return ret;
};

//----------------------------------------------------------------------
// Divide a tuple by a scalar
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S>
std::tuple<T1, T2, T3> operator/(const std::tuple<T1, T2, T3>& a, const S b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) / b,
			      std::get<1>(a) / b,
			      std::get<2>(a) / b );
  return ret;
};

//----------------------------------------------------------------------
// Add two tuples (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator+=(std::tuple<T1, T2, T3>& a, const std::tuple<T1, T2, T3>& b)
{
  std::get<0>(a)  += std::get<0>(b);
  std::get<1>(a)  += std::get<1>(b);
  std::get<2>(a)  += std::get<2>(b);

  return a;
};

//----------------------------------------------------------------------
// Multiple a tuple by a scalar (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S >
std::tuple<T1, T2, T3> operator*=(std::tuple<T1, T2, T3>& a, const S b)
{
  std::get<0>(a) *= b;
  std::get<1>(a) *= b;
  std::get<2>(a) *= b;

  return a;
};

//----------------------------------------------------------------------
// Divide a tuple by a scalar (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S >
std::tuple<T1, T2, T3> operator/=(std::tuple<T1, T2, T3>& a, const S b)
{
  std::get<0>(a) /= b;
  std::get<1>(a) /= b;
  std::get<2>(a) /= b;

  return a;
};

//----------------------------------------------------------------------
// Subtract two tuples (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator-=(std::tuple<T1, T2, T3>& a, const std::tuple<T1, T2, T3>& b)
{
  std::get<0>(a) -= std::get<0>(b);
  std::get<1>(a) -= std::get<1>(b);
  std::get<2>(a) -= std::get<2>(b);

  return a;
};

//----------------------------------------------------------------------
// Subtract two tuples (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::ostream& operator<< (std::ostream& stream, 
			  const std::tuple<T1, T2, T3>& v)
{
  stream << "["  << std::get<0>(v) 
	 << ", " << std::get<1>(v) 
	 << ", " << std::get<2>(v) << "]";

  return stream;
}

//----------------------------------------------------------------------
// Cell class definition follows
//----------------------------------------------------------------------
//
class sCell
{
public:
  typedef std::tuple<double, double, double> vcTup;

protected:
  typedef std::deque<vcTup> dqTup;

  std::map<sKeyPair, dqTup>     VelCrsList;
  std::map<sKeyPair, vcTup>     VelCrsSum;
  std::map<sKeyPair, unsigned>  VelCrsNum;
  
public:
  
  // Maximum number of cached relative velocities
  static size_t VelCrsSZ; 
  
  static double VelCrsMin;

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
  vcTup VelCrsAvg(sKeyPair indx);
  std::map<sKeyPair, vcTup> VelCrsAvg();
  //@}
  
  //@{
  //! Add maximum to list
  void VelCrsAdd(std::map<sKeyPair, vcTup>& vals);
  void VelCrsAdd(sKeyPair indx, const vcTup& val);
  void VelCrsAdd(sKeyPair indx, double vrm, double crm) 
  { VelCrsAdd(indx, vcTup(vrm*crm, vrm, crm)); }
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
