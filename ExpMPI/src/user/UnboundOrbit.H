// This may look like C code, but it is really -*- C++ -*-

#ifndef _UnboundOrbit_H
#define _UnboundOrbit_h

#include <string>

#include <Vector.h>
#include <massmodel.h>
#include <model3d.h>

#include <ParamDatabase.H>
#include <Trajectory.H>

//! Computes an satellite orbits and tidal forces in fixed halo
class UnboundOrbit : public Trajectory
{
private:

  SphericalModelTable *m;
  AxiSymModel *model;

  vector<double> R, T, PHI;
  vector<double> Time, Xpos, Ypos, Zpos;
  double TOFFSET;

  Matrix rotate, rotateI;

				// Private members

  void parse_args(void);	// Set parameters from parmFile
  void set_parm(string& word, string& alu);

public:

  //! Constructor
  UnboundOrbit(const string &file);

  //! Destructor
  ~UnboundOrbit(void);

				// Members

  //! Get satellite position in halo frame
  Vector get_satellite_orbit(double T);

  //! Get satellite position in halo frame
  void get_satellite_orbit(double T, double *v);

};

#endif
