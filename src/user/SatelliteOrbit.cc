/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Follow satellite orbit
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
 *  major rewrite: incorpated in to SatelliteOrbit class 11/15/98
 *
 ****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <cmath>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <orbit.H>
#include <massmodel.H>

#include <model2d.H>
#include <model3d.H>
#include <isothermal.H>
#include <hernquist.H>
#include <mestel.H>
#include <toomre.H>
#include <exponential.H>

#include <localmpi.H>
#include <global.H>
#include <SatelliteOrbit.H>

#include <YamlCheck.H>
#include <EXPException.H>
				// External prototype for Euler matrix

Eigen::Matrix3d return_euler(double PHI, double THETA, double PSI, int BODY);

const std::set<std::string>
SatelliteOrbit::valid_keys = {
  "HALO_MODEL",
  "PERI", 
  "APO",       
  "RSAT",
  "INCLINE",
  "PSI",    
  "PHIP",
  "VROT",
  "RCORE",
  "RMODMIN",
  "RMODMAX",
  "RA",    
  "DIVERGE",
  "MAXIT",    
  "NUMDF",
  "DIVRG_RFAC",
  "MODFILE", 
  "CIRCULAR",
  "orbfile",   
  "orbtmin",
  "orbtmax",
  "orbdelt",
};

// ===================================================================
// Constructor
// ===================================================================

SatelliteOrbit::SatelliteOrbit(const YAML::Node& conf)
{
  // Default parameters values
  //
  int    HALO_MODEL   = 0;
  double PERI         = 0.5;
  double APO          = 1.0;
  double RSAT         = 1.0;
  double INCLINE      = 45.0;
  double PSI          = 0.0;
  double PHIP         = 0.0;
  double VROT         = 1.0;
  double RCORE        = 1.0;
  double RMODMIN      = 1.0e-3;
  double RMODMAX      = 20.0;
  double RA           = 1.0e20;
  int    DIVERGE      = 0;
  int    MAXIT        = 2000;
  int    NUMDF        = 800;
  double DIVRG_RFAC   = 1.0;
  bool   orbfile      = true;
  double orbtmin      = -2.0;
  double orbtmax      =  2.0;
  double orbdelt      =  0.1;
  /*  */ circ         = false;
  std::string
         MODFILE      = "halo.model";


  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("SatelliteOrbit", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Get configured values
  //
  try {
    if (conf["HALO_MODEL"])  HALO_MODEL   = conf["HALO_MODEL"].as<int>();
    if (conf["PERI"      ])  PERI         = conf["PERI"      ].as<double>();
    if (conf["APO"       ])  APO          = conf["APO"       ].as<double>();
    if (conf["RSAT"      ])  RSAT         = conf["RSAT"      ].as<double>();
    if (conf["INCLINE"   ])  INCLINE      = conf["INCLINE"   ].as<double>();
    if (conf["PSI"       ])  PSI          = conf["PSI"       ].as<double>();
    if (conf["PHIP"      ])  PHIP         = conf["PHIP"      ].as<double>();
    if (conf["VROT"      ])  VROT         = conf["VROT"      ].as<double>();
    if (conf["RCORE"     ])  VROT         = conf["RCORE"     ].as<double>();
    if (conf["RMODMIN"   ])  RMODMIN      = conf["RMODMIN"   ].as<double>();
    if (conf["RMODMAX"   ])  RMODMAX      = conf["RMODMAX"   ].as<double>();
    if (conf["RA"        ])  RA           = conf["RA"        ].as<double>();
    if (conf["DIVERGE"   ])  DIVERGE      = conf["DIVERGE"   ].as<int>();
    if (conf["MAXIT"     ])  MAXIT        = conf["MAXIT"     ].as<int>();
    if (conf["NUMDF"     ])  NUMDF        = conf["NUMDF"     ].as<int>();
    if (conf["DIVRG_RFAC"])  DIVRG_RFAC   = conf["DIVRG_RFAC"].as<double>();
    if (conf["MODFILE"   ])  MODFILE      = conf["MODFILE"   ].as<std::string>();
    if (conf["CIRCULAR"  ])  circ         = conf["CIRCULAR"  ].as<bool>();
    if (conf["orbfile"   ])  orbfile      = conf["orbfile"   ].as<bool>();
    if (conf["orbtmin"   ])  orbtmin      = conf["orbtmin"   ].as<double>();
    if (conf["orbtmax"   ])  orbtmax      = conf["orbtmax"   ].as<double>();
    if (conf["orbdelt"   ])  orbdelt      = conf["orbdelt"   ].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in SatelliteOrbit: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

				// Initilize HALO model
  switch (HALO_MODEL) {
  case file:
    m = std::make_shared<SphericalModelTable>(MODFILE, DIVERGE, DIVRG_RFAC);
    m->setup_df(NUMDF, RA);
    halo_model = m;
				// Assign filename to ID string
    Model3dNames[0] = MODFILE;
    break;
    
  case sing_isothermal:
    halo_model = std::make_shared<SingIsothermalSphere>(VROT, RMODMIN, RMODMAX);
    break;

  case isothermal:
    halo_model = std::make_shared<IsothermalSphere>(RCORE, RMODMAX, VROT);
    break;

  case hernquist_model:
    halo_model = std::make_shared<HernquistSphere>(1.0, RMODMIN, RMODMAX);
    break; 

  default:
    cerr << "Illegal HALO model: " << HALO_MODEL << '\n';
    exit(-1);
  }
  

  // ===================================================================
  // Setup orbit
  // ===================================================================
    
  INCLINE *= M_PI/180.0;
  PSI     *= M_PI/180.0;
  PHIP    *= M_PI/180.0;
  
				// Set orientation of satellite orbit
  rotate  = return_euler(PHIP, INCLINE, PSI, 1);
  rotateI = rotate.inverse();
  
				// In case non-inertial is not desired
  omega = domega = 0.0;

  if (circ) {

    rsat = RSAT;
    vsat = sqrt(rsat*halo_model->get_dpot(rsat));
    Omega = vsat/rsat;

  } else {

    FindOrb::MAXIT = MAXIT;

    orb = std::make_shared<FindOrb>(halo_model, PERI, APO);

    OrbValues ret = orb->Anneal();

    if (myid==0) {
      cout << left << setw(60) << setfill('-') << '-' << endl << setfill(' ')
	   << setw(30) << "Boltzmann constant" 
	   << " | " << ret.Boltzmann << endl
	   << setw(30) << "Initial temperature" 
	   << " | " << ret.t0 << endl
	   << setw(30) << "Final temperature" 
	   << " | " << ret.tf << endl
	   << setw(30) << "Estimated minumum" 
	   << " | " << ret.energy << ", " << ret.kappa << endl
	   << setw(30) << "Functional value" 
	   << " | " << ret.value << endl
	   << setw(30) << "Peri, apo" 
	   << " | " << ret.peri << ", " << ret.apo << endl
	   << setw(30) << "Radial period" 
	   << " | " << ret.radial_period << endl
	   << setw(30) << "Azimuthal period" 
	   << " | " << ret.azimuthal_period << endl;
    }
  }
    
  if (myid==0 and orbfile) {	// Print diagnostic orbit file

    std::ofstream out(runtag + ".xyz");

    if (out.good()) {

      double T = orbtmin;

      while (T < orbtmax) {
	auto ps = get_satellite_orbit(T);

	out << setw(18) << T
	    << setw(18) << ps[0]
	    << setw(18) << ps[1]
	    << setw(18) << ps[2]
	    << std::endl;
	
	T += orbdelt;
      }
    }
    // END: orbit output
  }

}

// ===================================================================
// Destructor
// ===================================================================

SatelliteOrbit::~SatelliteOrbit(void)
{
  // Nothing
}

Eigen::Vector3d SatelliteOrbit::get_satellite_orbit(double T)
{
  double r, phi;

  if (circ) {
    r = rsat;
    phi = Omega*T;
  } else {
    r = orb->Orb().get_angle(6, T);
    phi = orb->Orb().get_angle(7, T);
  }

  v0[0] = r*cos(phi);
  v0[1] = r*sin(phi);
  v0[2] = 0.0;

				// Set current satellite position
  currentTime = T;
  currentR = rotate*v0;

  return currentR;
}

void SatelliteOrbit::get_satellite_orbit(double T, double *v)
{
  double r, phi;

  if (circ) {
    r = rsat;
    phi = Omega*T;
  } else {
    r = orb->Orb().get_angle(6, T);
    phi = orb->Orb().get_angle(7, T);
  }

  v0[0] = r*cos(phi);
  v0[1] = r*sin(phi);
  v0[2] = 0.0;

				// Set current satellite position
  currentTime = T;
  currentR = rotate*v0;

  for (int k=0; k<3; k++) v[k] = currentR[k];
}

