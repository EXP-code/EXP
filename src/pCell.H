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


//----------------------------------------------------------------------
// Cell class definition follows
//----------------------------------------------------------------------
//
class sCell
{
public:

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
  map<speciesKey,vector<double> > state;
  vector<double> stotal;

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
  std::set<unsigned long> bodies;
  
  //! For debugging . . . 
  void match(pCell* target, int& mcount, key_type& key);

  //! Collect bodies
  void collect(std::set<unsigned long>& bset);

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
  
  //! Attributes
  std::map<std::string, double> dattrib;
  std::map<std::string, int>    iattrib;

  key_type mykey;
  key_type mask;
  unsigned maxplev;
  
  pCell* parent;
  pCell* sample;
  map<unsigned, pCell*> children;
  bool isLeaf;
  
  key_indx keys;
  std::vector<unsigned long> bods;
  
  //! Track last time cell was evaluated
  double time;

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
  
  //! Find the cells whose body number is "just" larger than Bucket.
  //! Optional "message" string is for tagging debugging lines
  pCell* findSampleCell(const std::string& message="");
  
  //@{
  //! Return bodies from this cell
  Particle* Body(vector<unsigned long>::iterator k);
  Particle* Body(unsigned long i);
  //@}
  
  //! Return velocity statistics for particles in the cell
  void Vel(double &mass, vector<double>& v, vector<double>& v2);
  
  //! Return total volume of cell
  double Volume();
  
  //! Return the length scale factor from the quad tree
  double Scale();
  
  //! Return list of bodies in the leaves with this cell as a parent
  std::set<unsigned long>& Bodies()
  {
    bodies.clear();
    collect(bodies);
    return bodies;
  }

  //! For debugging: test that this cell is in the sample cell
  bool sampleTest()
  {
    int mcount = 0;
    key_type key = -1;
    sample->match(this, mcount, key);
    return mcount == 1;
  }

};

std::ostream& operator<< (std::ostream& stream, const sKeyPair& v);

#endif
