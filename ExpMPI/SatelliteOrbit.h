// This may look like C code, but it is really -*- C++ -*-

// ======================================================================
// Class to compute orbit of a satellite in spherical halo
//
//	o Orientation of orbit specified by Euler angles
//
// 	o Returns position and force in halo frame
//
//      o Compute tidal force in satellite frame 
//
//	o Arbitrary orientation of satellite body specified by Euler angles
//
// ======================================================================

#ifndef _SatelliteOrbit_h
#ifdef __GNUG__
#pragma interface
#endif
#define _SatelliteOrbit_h

#include <string>

#include <vector.h>
#include <massmodel.h>
#include <model3d.h>

class SatelliteOrbit
{
private:

  SphericalModelTable *m;
  AxiSymModel *halo_model;
  SphericalOrbit orb;
  Matrix rotate, rotateI;
  Three_Vector v, v0, u, u0, non;
  Matrix tidalRot, tidalRotI;
  double omega, domega;

  int halo_type;

				// Keep current satellte position
  double currentTime;
  Three_Vector currentR, currentF;

				// Private members

  void parse_args(void);	// Set parameters from parmFile
  void set_parm(char *word, char *valu);

				// General function for double and Vector input
  Vector get_tidal_force();
  Vector get_tidal_force_non_inertial();

				// Default initial parameters
  static Models3d HALO_MODEL;
  static double E;
  static double K;
  static double INCLINE;
  static double PSI;
  static double PHIP;
  static double VROT;
  static double RCORE;
  static double RMODMIN;
  static double RMODMAX;
  static double RA;
  static int DIVERGE;
  static int NUMDF;
  static int NRECS;
  static double DIVERGE_RFAC;
  static String MODFILE;

public:
  static string paramFile;

				// Constructor (no arguments);
  SatelliteOrbit(void);
				// Destructor
  ~SatelliteOrbit(void);

				// Members

				// Get satellite position in halo frame
  Vector get_satellite_orbit(double T);
				// Get force on satellite in halo frame
  Vector get_satellite_force(double T);


				// Member functions
				// for TIDAL calculation

				// Call once to set satelliet body orientation
  void setTidalOrientation(double phi, double theta, double psi);
				// Call to set satellite position
  void setTidalPosition(double T, int NI=0);

				// Retrieve satellite time
  double Time(void) { return currentTime; }

				// Get tidal force
  Vector tidalForce(const Vector p);
  Vector tidalForce(const double x, const double y, const double z);
  Vector tidalForce(const Vector p, const Vector q);
  Vector tidalForce(const double x, const double y, const double z,
		    const double u, const double v, const double w);

};

#endif
