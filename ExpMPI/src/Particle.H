
#ifndef Particle_H
#define Particle_H

#include <vector>
#include <map>
#ifdef I128
#include <uint128.H>
#endif

using namespace std;

//! Keeps track of all info for one particle
/*!
  The iattrib and dattrib vectors are used by individual components to
  carry additional parameters specific to different particle types.
 */
class Particle
{
public:

  //! Default effort value
  static float effort_default;

  //! Mass of particle
  double mass;

  //! Position vector
  double pos[3];

  //! Velocity vector
  double vel[3];

  //! Acceleration vector
  double acc[3];

  //! Gravitational potential due to particles  
  double pot;

  //! Gravitational potential from external forces
  double potext;
  
  //! Integer attributes
  vector<int> iattrib;

  //! Real (double) attributes
  vector<double> dattrib;

  //! Multistep level
  unsigned level;

  //! Internal process time step (<0 means ignore)
  float dtreq;

  //! Internal size scale (<0 means ignore)
  float scale;

  //! Scale factor for partitioning work
  float effort;

  //! Index for ordering
  unsigned long indx;

  //! Tree key
  unsigned tree;

  //! Hash key
#ifdef I128
  uint128 key;
#else
  unsigned long key;
#endif

  //! Constructor
  Particle();

  //! Constructor with presized attribute lists
  Particle(unsigned niatr, unsigned ndatr);

  //! Copy constructor
  Particle(const Particle &);

};

typedef map<unsigned long, Particle> PartMap;
typedef map<unsigned long, Particle>::iterator PartMapItr;


#endif
