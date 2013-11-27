// This may look like C code, but it is really -*- C++ -*-

// ======================================================================
// Class to compute orbit of a satellite in spherical halo
//
//	o Orientation of orbit specified by Euler angles
//
// 	o Returns position and force in halo frame
//
//	o Arbitrary orientation of satellite body specified by Euler angles
//
// ======================================================================

#ifndef _SatelliteOrbit_h
#define _SatelliteOrbit_h

#include <string>

#include <Vector.h>
#include <massmodel.h>
#include <model3d.h>

#include <FindOrb.H>
#include <ParamDatabase.H>
#include <Trajectory.H>

//! Computes an satellite orbits and tidal forces in fixed halo
class SatelliteOrbit : public Trajectory
{
private:

  SphericalModelTable *m;
  AxiSymModel *halo_model;
  FindOrb *orb;
  Matrix rotate, rotateI;
  Three_Vector v, v0, u, u0, non;
  Matrix tidalRot, tidalRotI;
  double omega, domega;
  double rsat, vsat, Omega;
  bool circ;

  int halo_type;

				// Keep current satellte position
  double currentTime;
  Three_Vector currentR;

public:

  //! Constructor
  SatelliteOrbit(const string &file);

  //! Destructor
  ~SatelliteOrbit(void);

				// Members

  //! Get satellite position in halo frame
  Vector get_satellite_orbit(double T);

  //! Get satellite position in halo frame
  void get_satellite_orbit(double T, double *v);

};

#endif
