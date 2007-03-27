
#ifndef Particle_H
#define Particle_H

#include <vector>

using namespace std;

//! Keeps track of all info for one particle
/*!
  The iattrib and dattrib vectors are used by individual components to
  carry additional parameters specific to different particle types.
 */
class Particle
{
public:

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
  vector<int>    iattrib;

  //! Real (double) attributes
  vector<double> dattrib;

  //! Index for ordering
  unsigned long indx;

  //! Constructor
  Particle();

  //! Copy constructor
  Particle(const Particle &);

};

#endif
