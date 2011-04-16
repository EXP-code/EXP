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
 *  major rewrite: incorpated in to SatelliteOrbit class 11/15/98
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

  Rperi = config->get<double>("Rperi");
  Vsat  = config->get<double>("Vsat");

  double THETA   = config->get<double>("THETA")   * M_PI/180.0;
  double PSI     = config->get<double>("PSI")     * M_PI/180.0;
  double PHIP    = config->get<double>("PHIP")    * M_PI/180.0;

  rotate  = return_euler_slater(PHIP, THETA, PSI, 1);

  if (myid==0) {
    
    Vector ret(1, 3);
    ret[1] = config->get<double>("X0"); 
    ret[2] = config->get<double>("Y0");
    ret[3] = config->get<double>("Z0"); 

    cout << "LinearOrbit initiated with:" << endl
	 << setw(10) << "" << setw(10) << "THETA" 
	 << " = " << config->get<double>("THETA") << endl
	 << setw(10) << "" << setw(10) << "PSI" 
	 << " = " << config->get<double>("PSI") << endl
	 << setw(10) << "" << setw(10) << "PHIP"
	 << " = " << config->get<double>("PHIP") << endl
	 << "Initial position and velocity is:" << endl
	 << setw(10) << "" << setw(10) << "(X, Y, Z)" 
	 << " = (" << ret[1] 
	 << ", "   << ret[2]
	 << ", "   << ret[3]
	 << ")" << endl
	 << setw(10) << "" << setw(10) << "(U, V, W)" 
	 << " = (" << 0
	 << ", "   << Vsat
	 << ", "   << 0
	 << ")" << endl
	 << "Rotated position and velocity is:" << endl;

    ret = rotate * ret;
    cout << setw(10) << "" << setw(10) << "(X, Y, Z)" 
	 << " = (" << ret[1] 
	 << ", "   << ret[2]
	 << ", "   << ret[3]
	 << ")" << endl;

    ret[1] = ret[3] = 0.0;
    ret[2] = Vsat;
    ret = rotate * ret;
    cout << setw(10) << "" << setw(10) << "(U, V, W)" 
	 << " = (" << ret[1]
	 << ", "   << ret[2]
	 << ", "   << ret[3]
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

  ret[1] = Rperi;
  ret[2] = Vsat * t;
  ret[3] = 0.0;

  ret = rotate * ret;

  return ret;
}

void LinearOrbit::get_satellite_orbit(double t, double *v)
{
  Vector ret(1, 3);

  ret[1] = 0.0;
  ret[2] = Vsat; 
  ret[3] = 0.0;

  ret = rotate * ret;
  v[0] = ret[1];
  v[1] = ret[2];
  v[2] = ret[3];
}

