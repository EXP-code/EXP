// This may look like C code, but it is really -*- C++ -*-

#ifndef _LinearOrbit_H
#define _LinearOrbit_H

#include <string>

#include <Trajectory.H>

//! Computes a linear satellite trajectory
class LinearOrbit : public Trajectory
{
private:

  double X0, Y0, Z0, Vsat;
  Eigen::Matrix3d rotate;

				// Private members

  void parse_args(void);	// Set parameters from parmFile
  void set_parm(string& word, string& alu);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  LinearOrbit(const YAML::Node& conf);

  //! Destructor
  ~LinearOrbit(void);

				// Members

  //! Get satellite position in halo frame
  Eigen::Vector3d get_satellite_orbit(double T);

  //! Get satellite position in halo frame
  void get_satellite_orbit(double T, double *v)
  {
    auto ret = get_satellite_orbit(T);
    for (unsigned i=0; i<3; i++) v[i] = ret[i];
  }

};

#endif
