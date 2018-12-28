#ifndef _Trajectory_H
#define _Trajectory_H

#include <string>

#include <Vector.h>
#include <massmodel.h>
#include <model3d.h>

#include <FindOrb.H>

//! Computes an satellite orbits and tidal forces in fixed halo
class Trajectory
{
protected:

public:

  //! Constructor
  Trajectory() {}

  //! Destructor
  virtual ~Trajectory() {}

				// Members
  //! Get satellite position
  virtual Vector get_satellite_orbit(double T) = 0;

  //! Get satellite position
  virtual void get_satellite_orbit(double T, double *v) = 0;
};

#endif
