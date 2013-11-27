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

#include <kevin_complex.h>
#include <Vector.h>
#include <orbit.h>
#include <LinearOrbit.H>

using namespace std;

#include <global.H>
				// External prototype for Euler matrix

Matrix return_euler_slater(double PHI, double THETA, double PSI, int BODY);

				// Default input parameters (storage)

static database_record init[] = {
  {"Vsat",	"double",	"1.0"},
  {"X0",	"double",	"0.0"},
  {"Y0",	"double",	"0.0"},
  {"Z0",	"double",	"0.0"},
  {"PHIP",	"double",	"0.0"},
  {"THETA",	"double",	"0.0"},
  {"PSI",	"double",	"0.0"},
  {""		"",		""}
};

// ===================================================================
// Constructor
// ===================================================================

LinearOrbit::LinearOrbit(const string &conf)
{
  config = new ParamDatabase(init);
  config->parseFile(conf);

// =================================
// Initialize parameters
// =================================

  X0 = config->get<double>("X0");
  Y0 = config->get<double>("Y0");
  Z0 = config->get<double>("Z0");
  Vsat  = config->get<double>("Vsat");

  double THETA   = config->get<double>("THETA")   * M_PI/180.0;
  double PSI     = config->get<double>("PSI")     * M_PI/180.0;
  double PHIP    = config->get<double>("PHIP")    * M_PI/180.0;

  rotate  = return_euler_slater(PHIP, THETA, PSI, 1);

  if (myid==0) {
    
    Vector a(1, 3), b(1, 3);
    a[1] = config->get<double>("X0"); 
    a[2] = config->get<double>("Y0");
    a[3] = config->get<double>("Z0"); 

    cout << "LinearOrbit initialized with:" << endl << left
	 << setw(5) << "" << setw(10) << "THETA" << " = "
	 << config->get<double>("THETA") << endl
	 << setw(5) << "" << setw(10) << "PSI "  << " = " 
	 << config->get<double>("PSI")   << endl
	 << setw(5) << "" << setw(10) << "PHIP"  << " = " 
	 << config->get<double>("PHIP") << endl
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
// Destructior
// ===================================================================

LinearOrbit::~LinearOrbit(void)
{
  delete config;
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

