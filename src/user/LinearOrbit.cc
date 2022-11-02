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

#include <YamlCheck.H>
#include <EXPException.H>
#include <orbit.H>
#include <LinearOrbit.H>

using namespace std;

#include <global.H>
				// External prototype for Euler matrix

Eigen::Matrix3d
return_euler_slater(double PHI, double THETA, double PSI, int BODY);

const std::set<std::string>
LinearOrbit::valid_keys = {
  "X0",
  "Y0",
  "Z0",
  "Vsat",
  "THETA",
  "PSI",
  "PHIP"
};


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

  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("LinearOrbit", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
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
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  rotate = return_euler_slater(PHIP, THETA, PSI, 1);

  if (myid==0) {
    
    Eigen::Vector3d a, b;
    a[0] = X0;
    a[1] = Y0;
    a[2] = Z0;

    cout << "LinearOrbit initialized with:" << endl << left
	 << setw(5) << "" << setw(10) << "THETA" << " = "
	 << THETA * 180.0/M_PI << endl
	 << setw(5) << "" << setw(10) << "PSI "  << " = " 
	 << PSI * 180.0/M_PI   << endl
	 << setw(5) << "" << setw(10) << "PHIP"  << " = " 
	 << PHIP * 180.0/M_PI << endl
	 << "Initial position and velocity is:" << endl
	 << setw(5) << "" << setw(10) << "(X, Y, Z)" 
	 << " = (" << a[0]
	 << ", "   << a[1]
	 << ", "   << a[2]
	 << ")" << endl
	 << setw(5) << "" << setw(10) << "(U, V, W)" 
	 << " = (" << 0
	 << ", "   << Vsat
	 << ", "   << 0
	 << ")" << endl
	 << "Rotated position and velocity is:" << endl;

    b = rotate * a;
    cout << setw(5) << "" << setw(10) << "(X, Y, Z)" 
	 << " = (" << b[0] 
	 << ", "   << b[1]
	 << ", "   << b[2]
	 << ")" << endl;

    a[0] = a[2] = 0.0;
    a[1] = Vsat;
    b = rotate * a;
    cout << setw(5) << "" << setw(10) << "(U, V, W)" 
	 << " = (" << b[0]
	 << ", "   << b[1]
	 << ", "   << b[2]
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

Eigen::Vector3d LinearOrbit::get_satellite_orbit(double t)
{
  Eigen::Vector3d ret;

  ret[0] = X0;
  ret[1] = Y0 + Vsat * t;
  ret[2] = Z0;

  ret = rotate * ret;

  return ret;
}

