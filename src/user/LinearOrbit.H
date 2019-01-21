// This may look like C code, but it is really -*- C++ -*-

#ifndef _LinearOrbit_H
#define _LinearOrbit_h

#include <string>

#include <Vector.h>
#include <Trajectory.H>

//! Computes a linear satellite trajectory
class LinearOrbit : public Trajectory
{
private:

  double X0, Y0, Z0, Vsat;
  Matrix rotate;

				// Private members

  void parse_args(void);	// Set parameters from parmFile
  void set_parm(string& word, string& alu);

public:

  //! Constructor
  LinearOrbit(const YAML::Node& conf);

  //! Destructor
  ~LinearOrbit(void);

				// Members

  //! Get satellite position in halo frame
  Vector get_satellite_orbit(double T);

  //! Get satellite position in halo frame
  void get_satellite_orbit(double T, double *v)
  {
    Vector ret = get_satellite_orbit(T);
    for (unsigned i=0; i<3; i++) v[i] = ret[i+1];
  }

};

#endif
