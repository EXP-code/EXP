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

#include <stdlib.h>
#include <math.h>
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

#include <SatelliteOrbit.h>

#include <localmpi.h>
#include "global.H"
				// External prototype for Euler matrix

Matrix return_euler(double PHI, double THETA, double PSI, int BODY);

				// Default input parameters (storage)

database_record init[] = {
  {"HALO_MODEL",	"int",		"0"},
  {"PERI",		"double",	"0.5"},
  {"APO",		"double",	"1.0"},
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
  {"NRECS",		"int", 		"512"},
  {"DIVERGE_RFAC",	"double",	"1.0"},
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
				config->get<double>("DIVERGE_RFAC"));
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

  double INCLINE = config->get<double>("INCLINE");
  double PSI     = config->get<double>("PSI");
  double PHIP    = config->get<double>("PHIP");

  INCLINE *= M_PI/180.0;
  PSI     *= M_PI/180.0;
  PHIP    *= M_PI/180.0;
  
				// Set orientation of satellite orbit
  rotate  = return_euler(PHIP, INCLINE, PSI, 1);
  rotateI = rotate.Inverse();

				// Set default body rotation to identity
  setTidalOrientation(0.0, 0.0, 0.0);

				// In case non-inertial is not desired
  omega = domega = 0.0;

  if (myid==0) {

    double r   = orb->Orb().get_angle(6, 0.0);
    double phi = orb->Orb().get_angle(7, 0.0);

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
				// Delete halo model (switch to
				// call correct desstructor)
  switch (halo_type) {
  case file:
    delete m;
    break;
    
  case sing_isothermal:
    delete (SingIsothermalSphere *)halo_model;
    break;

  case isothermal:
    delete (IsothermalSphere *)halo_model;
    break;

  case hernquist_model:
    delete (HernquistSphere *)hernquist_model;
    break; 
  }
  
  delete orb;

}

// ===================================================================
// Set satellite body transformations
// ===================================================================

void SatelliteOrbit::setTidalOrientation(double phi, double theta, double psi)
{
				// Transformation from body to halo coordinates
  tidalRot = return_euler(phi, theta, psi, 1);
				// Transformation back to body coordinates
  tidalRotI = tidalRot.Inverse();
}


// ===================================================================
// The first two package position vector for tidal force computation
// by get_tidal_force()
// ===================================================================

Vector SatelliteOrbit::tidalForce(const Vector p)
{
  v0[1] = p[1];
  v0[2] = p[2];
  v0[3] = p[3];

  return get_tidal_force();
}

Vector SatelliteOrbit::tidalForce(const double x, 
				  const double y,
				  const double z)
{
  v0[1] = x;
  v0[2] = y;
  v0[3] = z;

  return get_tidal_force();
}

Vector SatelliteOrbit::tidalForce(const Vector p, const Vector q)
{
  v0[1] = p[1];
  v0[2] = p[2];
  v0[3] = p[3];

  u0[1] = q[1];
  u0[2] = q[2];
  u0[3] = q[3];

  return get_tidal_force_non_inertial();
}

Vector SatelliteOrbit::tidalForce(const double x, 
				  const double y,
				  const double z,
				  const double u,
				  const double v,
				  const double w)
{
  v0[1] = x;
  v0[2] = y;
  v0[3] = z;

  u0[1] = u;
  u0[2] = v;
  u0[3] = w;

  return get_tidal_force_non_inertial();
}

Vector SatelliteOrbit::get_tidal_force(void)
{
  v = currentR + tidalRot * v0;

  double r = sqrt(v*v);
  double dpot = halo_model->get_dpot(r);
  
  v0 = -dpot*v/r - currentF;

  return tidalRotI * v0;
}


Vector SatelliteOrbit::get_tidal_force_non_inertial(void)
{
  v = rotateI * tidalRot * v0;
  u = rotateI * tidalRot * u0;
  //  ^         ^
  //  |         |
  //  |         \-- go to halo coordinates
  //  |           
  //  \------------ go from halo to orbital plane coordinates
  //

				// Non inertial forces in orbital plane

  non[1] =  2.0*omega*u[2] + omega*omega*v[1] + domega*v[2];
  non[2] = -2.0*omega*u[1] + omega*omega*v[2] - domega*v[1];
  non[3] = 0.0;

				// Compute force vector in halo coordinates

  v = currentR + tidalRot * v0;
  //             ^
  //             |
  //             \-- go to halo coordinates


  double r = sqrt(v*v);
  double dpot = halo_model->get_dpot(r);
  
  v0 = -dpot*v/r - currentF + rotate*non;
  //                          ^
  //                          |           
  //                          \--- go from orbital plane to halo coordinates
  //
  //
  //     /-- go from halo back to body coordinates
  //     |
  //     v
  return tidalRotI * v0;
}


void SatelliteOrbit::setTidalPosition(double T, int NI)
{
  double r = orb->Orb().get_angle(6, T);
  double phi = orb->Orb().get_angle(7, T);
  double dpot = halo_model->get_dpot(r);

  currentTime = T;

  v0[1] = r*cos(phi);
  v0[2] = r*sin(phi);
  v0[3] = 0.0;

				// Set current satellite position
  currentR = rotate*v0;

  v0[1] = -dpot*cos(phi);
  v0[2] = -dpot*sin(phi);
  v0[3] = 0.0;

				// Set com force
  currentF = rotate*v0;


  if (NI) {			// Set up for non-inertial terms

    double delT = 2.0*M_PI/orb->Orb().get_freq(2)/40.0;
  
    omega = (
	     orb->Orb().get_angle(7, T+0.5*delT) -
	     orb->Orb().get_angle(7, T-0.5*delT)
	     ) / delT;

    domega = (
	      orb->Orb().get_angle(7, T+delT) -
	      orb->Orb().get_angle(7, T     ) * 2.0 +
	      orb->Orb().get_angle(7, T-delT)
	      ) / (delT*delT);
  }
}


Vector SatelliteOrbit::get_satellite_orbit(double T)
{
  double r = orb->Orb().get_angle(6, T);
  double phi = orb->Orb().get_angle(7, T);

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
  double r = orb->Orb().get_angle(6, T);
  double phi = orb->Orb().get_angle(7, T);

  v0[1] = r*cos(phi);
  v0[2] = r*sin(phi);
  v0[3] = 0.0;

				// Set current satellite position
  currentTime = T;
  currentR = rotate*v0;

  for (int k=0; k<3; k++) v[k] = currentR[k+1];
}

Vector SatelliteOrbit::get_satellite_force(double T)
{
  double r = orb->Orb().get_angle(6, T);
  double phi = orb->Orb().get_angle(7, T);
  double dpot = halo_model->get_dpot(r);

  v0[1] = -dpot*cos(phi);
  v0[2] = -dpot*sin(phi);
  v0[3] = 0.0;

  return rotate*v0;
}
