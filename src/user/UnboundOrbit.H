// This may look like C code, but it is really -*- C++ -*-

#ifndef _UnboundOrbit_H
#define _UnboundOrbit_h

#include <memory>
#include <string>

#include <yaml-cpp/yaml.h>

#include <massmodel.H>
#include <model3d.H>

#include <Trajectory.H>

//! Computes an satellite orbits and tidal forces in fixed halo
class UnboundOrbit : public Trajectory
{
private:

  std::shared_ptr<SphericalModelTable> m;
  std::shared_ptr<AxiSymModel> model;

  vector<double> R, T, PHI;
  vector<double> Time, Xpos, Ypos, Zpos;
  double TOFFSET;

  Eigen::Matrix3d rotate, rotateI;

				// Private members

  void parse_args(void);	// Set parameters from parmFile
  void set_parm(string& word, string& alu);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UnboundOrbit(const YAML::Node& conf);

  //! Destructor
  ~UnboundOrbit(void);

				// Members

  //! Get satellite position in halo frame
  Eigen::Vector3d get_satellite_orbit(double T);

  //! Get satellite position in halo frame
  void get_satellite_orbit(double T, double *v);

};

#endif
