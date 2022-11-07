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

#ifndef _SatelliteOrbit_H
#define _SatelliteOrbit_H

#include <string>
#include <massmodel.H>
#include <model3d.H>

#include <FindOrb.H>
#include <Trajectory.H>

//! Computes an satellite orbits and tidal forces in fixed halo
class SatelliteOrbit : public Trajectory
{
private:

  std::shared_ptr<SphericalModelTable> m;
  std::shared_ptr<AxiSymModel> halo_model;
  std::shared_ptr<FindOrb> orb;
  Eigen::Matrix3d rotate, rotateI;
  Eigen::Vector3d v, v0, u, u0, non;
  Eigen::Matrix3d tidalRot, tidalRotI;
  double omega, domega;
  double rsat, vsat, Omega;
  bool circ;

  int halo_type;

				// Keep current satellte position
  double currentTime;
  Eigen::Vector3d currentR;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  SatelliteOrbit(const YAML::Node& conf);

  //! Destructor
  ~SatelliteOrbit(void);

				// Members

  //! Get satellite position in halo frame
  Eigen::Vector3d get_satellite_orbit(double T);

  //! Get satellite position in halo frame
  void get_satellite_orbit(double T, double *v);

};

#endif
