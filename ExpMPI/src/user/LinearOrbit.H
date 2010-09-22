// This may look like C code, but it is really -*- C++ -*-

#ifndef _LinearOrbit_H
#define _LinearOrbit_h

#include <string>

#include <Vector.h>
#include <ParamDatabase.H>
#include <Trajectory.H>

//! Computes a linear satellite trajectory
class LinearOrbit : public Trajectory
{
private:

  double Rperi, Vsat;
  Matrix rotate;

				// Private members

  void parse_args(void);	// Set parameters from parmFile
  void set_parm(string& word, string& alu);

public:

  //! Constructor
  LinearOrbit(const string &file);

  //! Destructor
  ~LinearOrbit(void);

				// Members

  //! Get satellite position in halo frame
  Vector get_satellite_orbit(double T);

  //! Get satellite position in halo frame
  void get_satellite_orbit(double T, double *v);

};

#endif
