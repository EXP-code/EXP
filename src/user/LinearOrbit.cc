// This may look like C code, but it is really -*- C++ -*-

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Follow linear satellite trajectory
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  revision MDW 01/10/93
 *  updated to use orbit classes
 *           MDW 07/15/94
 *
 *  major rewrite: incorporated in to SatelliteOrbit class 11/15/98
 *  changed position spec from impact parameter to full 3d vector
 *
 ****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

#include <yaml-cpp/yaml.h>

#include <Vector.h>
#include <orbit.h>
#include <LinearOrbit.H>

using namespace std;

#include <global.H>
				// External prototype for Euler matrix

Matrix return_euler_slater(double PHI, double THETA, double PSI, int BODY);

// ===================================================================
// Constructor
// ===================================================================

LinearOrbit::LinearOrbit(const YAML::Node& conf)
{
// =================================
// Initialize parameters
// =================================

  X0    = 0.0;
  Y0    = 0.0;
  Z0    = 0.0;
  Vsat  = 1.0;

  double PHIP  = 0.0;
  double THETA = 0.0;
  double PSI   = 0.0;

  try {
    if (conf["X0"])      X0      = conf["X0"   ].as<double>();
    if (conf["Y0"])      Y0      = conf["Y0"   ].as<double>();
    if (conf["Z0"])      Z0      = conf["Z0"   ].as<double>();
    if (conf["Vsat"])    Vsat    = conf["Vsat" ].as<double>();
    if (conf["THETA"])   THETA   = conf["THETA"].as<double>() * M_PI/180.0;
    if (conf["PSI"])     PSI     = conf["PSI"  ].as<double>() * M_PI/180.0;
    if (conf["PHIP"])    PHIP    = conf["PHIP" ].as<double>() * M_PI/180.0;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in LinearOrbit: "
			   << error.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  rotate = return_euler_slater(PHIP, THETA, PSI, 1);

  if (myid==0) {
    
    Vector a(1, 3), b(1, 3);
    a[1] = X0;
    a[2] = Y0;
    a[3] = Z0;

    cout << "LinearOrbit initialized with:" << endl << left
	 << setw(5) << "" << setw(10) << "THETA" << " = "
	 << THETA * 180.0/M_PI << endl
	 << setw(5) << "" << setw(10) << "PSI "  << " = " 
	 << PSI * 180.0/M_PI   << endl
	 << setw(5) << "" << setw(10) << "PHIP"  << " = " 
	 << PHIP * 180.0/M_PI << endl
	 << "Initial position and velocity is:" << endl
	 << setw(5) << "" << setw(10) << "(X, Y, Z)" 
	 << " = (" << a[1]
	 << ", "   << a[2]
	 << ", "   << a[3]
	 << ")" << endl
	 << setw(5) << "" << setw(10) << "(U, V, W)" 
	 << " = (" << 0
	 << ", "   << Vsat
	 << ", "   << 0
	 << ")" << endl
	 << "Rotated position and velocity is:" << endl;

    b = rotate * a;
    cout << setw(5) << "" << setw(10) << "(X, Y, Z)" 
	 << " = (" << b[1] 
	 << ", "   << b[2]
	 << ", "   << b[3]
	 << ")" << endl;

    a[1] = a[3] = 0.0;
    a[2] = Vsat;
    b = rotate * a;
    cout << setw(5) << "" << setw(10) << "(U, V, W)" 
	 << " = (" << b[1]
	 << ", "   << b[2]
	 << ", "   << b[3]
	 << ")" << endl;
  }

}


// ===================================================================
// Destructor
// ===================================================================

LinearOrbit::~LinearOrbit(void)
{
  // Nothing
}

Vector LinearOrbit::get_satellite_orbit(double t)
{
  Vector ret(1, 3);

  ret[1] = X0;
  ret[2] = Y0 + Vsat * t;
  ret[3] = Z0;

  ret = rotate * ret;

  return ret;
}

