// This may look like C code, but it is really -*- C++ -*-

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

// static char rcsid_SatelliteOrbit[] = "$Id$";

#include <cstdlib>
#include <cmath>
#include <string>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <kevin_complex.h>
#include <Vector.h>
#include <orbit.h>
#include <massmodel.h>

#include <model2d.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>
#include <mestel.h>
#include <toomre.h>
#include <exponential.h>

#include <localmpi.h>
#include <global.H>
#include <SatelliteOrbit.h>
#include <ParamDatabase.H>
				// External prototype for Euler matrix

Matrix return_euler(double PHI, double THETA, double PSI, int BODY);

				// Default input parameters (storage)

static database_record init[] = {
  {"HALO_MODEL",	"int",		"0"  },
  {"PERI",		"double",	"0.5"},
  {"APO",		"double",	"1.0"},
  {"RSAT",		"double",	"1.0"},
  {"INCLINE",		"double",	"45.0"},	  
  {"PSI",		"double",	"0.0"},
  {"PHIP",		"double",	"0.0"},
  {"VROT",		"double",	"1.0"},
  {"RCORE",		"double",	"1.0"},
  {"RMODMIN",		"double", 	"1.0e-3"},
  {"RMODMAX",		"double", 	"20.0"},
  {"RA",		"double", 	"1.0e20"},
  {"DIVERGE",		"int", 		"0"},
  {"NUMDF",		"int", 		"100"},
  {"MAXIT",		"int", 		"2000"},
  {"NRECS",		"int", 		"512"},
  {"DIVERGE_RFAC",	"double",	"1.0"},
  {"CIRCULAR",		"bool",		"false"},
  {"MODFILE",		"string",	"halo.model"},
  {"",			"",    		""}
};

// ===================================================================
// Constructor
// ===================================================================

SatelliteOrbit::SatelliteOrbit(const string &conf)
{
  config = new ParamDatabase(init);
  config->parseFile(conf);

				// Initilize HALO model
  switch (config->get<int>("HALO_MODEL")) {
  case file:
    m = new SphericalModelTable(
				config->get<string>("MODFILE"), 
				config->get<int>("DIVERGE"), 
				config->get<double>("DIVERGE_RFAC")
				);
    m->setup_df(config->get<int>("NUMDF"), config->get<double>("RA"));
    halo_model = m;
				// Assign filename to ID string
    Model3dNames[0] = config->get<string>("MODFILE");
    break;
    
  case sing_isothermal:
    halo_model = new SingIsothermalSphere(
					  config->get<double>("VROT"), 
					  config->get<double>("RMODMIN"),
					  config->get<double>("RMODMAX")
					  );
    break;

  case isothermal:
    halo_model = new IsothermalSphere(config->get<double>("RCORE"), 
				      config->get<double>("RMODMAX"), 
				      config->get<double>("VROT"));
    break;

  case hernquist_model:
    halo_model = new HernquistSphere(1.0, // Halo model
				     config->get<double>("RMODMIN"), 
				     config->get<double>("RMODMAX"));
    break; 

  default:
    cerr << "Illegal HALO model: " << config->get<int>("HALO_MODEL") << '\n';
    exit(-1);
  }
  

// ===================================================================
// Setup orbit
// ===================================================================

  double INCLINE = config->get<double>("INCLINE");
  double PSI     = config->get<double>("PSI");
  double PHIP    = config->get<double>("PHIP");
    
  INCLINE *= M_PI/180.0;
  PSI     *= M_PI/180.0;
  PHIP    *= M_PI/180.0;
  
				// Set orientation of satellite orbit
  rotate  = return_euler(PHIP, INCLINE, PSI, 1);
  rotateI = rotate.Inverse();
  
				// In case non-inertial is not desired
  omega = domega = 0.0;

  circ = config->get<bool>("CIRCULAR");

  if (circ) {

    rsat = config->get<double>("RSAT");
    vsat = sqrt(rsat*halo_model->get_dpot(rsat));
    Omega = vsat/rsat;

  } else {

    FindOrb::MAXIT = config->get<int>("MAXIT");

    orb = new FindOrb(
		      halo_model,
		      config->get<double>("PERI"), 
		      config->get<double>("APO")
		      );

    OrbValues ret = orb->Anneal();

    if (myid==0) {
      cout << left << setw(60) << setfill('-') << '-' << endl << setfill(' ');
      cout << "Boltzman constant: " << ret.Boltzmann << endl
	   << "Initial temperature: " << ret.t0 << '\t'
	   << "Final temperature: " << ret.tf << endl
	   << "Estimated minumum at: " << ret.energy
	   << ", " << ret.kappa << endl
	   << "Functional value = " << ret.value << endl
	   << "Peri, apo = " << ret.peri << ", " << ret.apo << endl
	   << "Radial period = " << ret.radial_period << endl
	   << "Aximuthal period = " << ret.azimuthal_period << endl;
    }
  }
    
  if (myid==0) {

    double r, phi;

    if (circ) {
      r = rsat;
      phi = 0.0;
    } else {
      r   = orb->Orb().get_angle(6, 0.0);
      phi = orb->Orb().get_angle(7, 0.0);
    }
      
    v0[1] = r*cos(phi);
    v0[2] = r*sin(phi);
    v0[3] = 0.0;
    
				// Set current satellite position
    currentR = rotate*v0;

    cout << "Position at T=0: x, y, z = " << currentR[1] << ", "
	 << currentR[2] << ", " << currentR[3] << endl
	 << setw(60) << setfill('-') << '-' << endl << setfill(' ');
  }

}

// ===================================================================
// Destructior
// ===================================================================

SatelliteOrbit::~SatelliteOrbit(void)
{
  if (m) delete m;
  else   delete halo_model;
  
  delete orb;
  delete config;
}

Vector SatelliteOrbit::get_satellite_orbit(double T)
{
  double r, phi;

  if (circ) {
    r = rsat;
    phi = Omega*T;
  } else {
    r = orb->Orb().get_angle(6, T);
    phi = orb->Orb().get_angle(7, T);
  }

  v0[1] = r*cos(phi);
  v0[2] = r*sin(phi);
  v0[3] = 0.0;

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

  v0[1] = r*cos(phi);
  v0[2] = r*sin(phi);
  v0[3] = 0.0;

				// Set current satellite position
  currentTime = T;
  currentR = rotate*v0;

  for (int k=0; k<3; k++) v[k] = currentR[k+1];
}

