// This may look like C code, but it is really -*- C++ -*-

#pragma implementation

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

#include "global.H"
				// External prototype for Euler matrix

Matrix return_euler(double PHI, double THETA, double PSI, int BODY);

				// Default input parameters (storage)

Models3d SatelliteOrbit::HALO_MODEL= (Models3d)file;
double SatelliteOrbit::E=-1.0;
double SatelliteOrbit::K=0.5;
double SatelliteOrbit::INCLINE=45.0;
double SatelliteOrbit::PSI=0.0;
double SatelliteOrbit::PHIP=0.0;
double SatelliteOrbit::VROT=1.0;
double SatelliteOrbit::RCORE=1.0;
double SatelliteOrbit::RMODMIN=1.0e-3;
double SatelliteOrbit::RMODMAX=20.0;
double SatelliteOrbit::RA=1.0e20;
int SatelliteOrbit::DIVERGE=0;
int SatelliteOrbit::NUMDF=100;
int SatelliteOrbit::NRECS=512;
double SatelliteOrbit::DIVERGE_RFAC=1.0;
string SatelliteOrbit::MODFILE="halo.model";
string SatelliteOrbit::paramFile = "satellite.params";


// ===================================================================
// Constructor
// ===================================================================

SatelliteOrbit::SatelliteOrbit(void)
{
  parse_args();

				// Initilize HALO model
  switch (HALO_MODEL) {
  case file:
    m = new SphericalModelTable(MODFILE, DIVERGE, DIVERGE_RFAC);
    m->setup_df(NUMDF, RA);
    halo_model = m;
    Model3dNames[0] = MODFILE;	// Assign filename to ID string
    break;
    
  case sing_isothermal:
    halo_model = new SingIsothermalSphere(VROT, RMODMIN, RMODMAX);
    break;

  case isothermal:
    halo_model = new IsothermalSphere(RCORE, RMODMAX, VROT);
    break;

  case hernquist_model:
    halo_model = new HernquistSphere(1.0, RMODMIN, RMODMAX); // Halo model
    break; 

  default:
    cerr << "Illegal HALO model: " << HALO_MODEL << '\n';
    exit(-1);
    }
  

// ===================================================================
// Setup orbit
// ===================================================================

  SphericalOrbit create(halo_model, E, K);
  create.set_numerical_params(NRECS);

  orb = create;

  INCLINE *= M_PI/180.0;
  PSI *= M_PI/180.0;
  PHIP *= M_PI/180.0;
  
				// Set orientation of satellite orbit
  rotate = return_euler(PHIP, INCLINE, PSI, 1);
  rotateI = rotate.Inverse();

				// Set default body rotation to identity
  setTidalOrientation(0.0, 0.0, 0.0);

				// In case non-inertial is not desired
  omega = domega = 0.0;

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
  double r = orb.get_angle(6, T);
  double phi = orb.get_angle(7, T);
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

    double delT = 2.0*M_PI/orb.get_freq(2)/40.0;
  
    omega = (
	     orb.get_angle(7, T+0.5*delT) -
	     orb.get_angle(7, T-0.5*delT)
	     ) / delT;

    domega = (
	      orb.get_angle(7, T+delT) -
	      orb.get_angle(7, T     ) * 2.0 +
	      orb.get_angle(7, T-delT)
	      ) / (delT*delT);
  }
}


Vector SatelliteOrbit::get_satellite_orbit(double T)
{
  double r = orb.get_angle(6, T);
  double phi = orb.get_angle(7, T);

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
  double r = orb.get_angle(6, T);
  double phi = orb.get_angle(7, T);

  v0[1] = r*cos(phi);
  v0[2] = r*sin(phi);
  v0[3] = 0.0;

				// Set current satellite position
  currentTime = T;
  currentR = rotate*v0;

  for (int k=0; k<3; k++) v[k] = currentR[k-1];
}

Vector SatelliteOrbit::get_satellite_force(double T)
{
  double r = orb.get_angle(6, T);
  double phi = orb.get_angle(7, T);
  double dpot = halo_model->get_dpot(r);

  v0[1] = -dpot*cos(phi);
  v0[2] = -dpot*sin(phi);
  v0[3] = 0.0;

  return rotate*v0;
}

// ===================================================================
// Initializatio helper functions
// ===================================================================


void SatelliteOrbit::set_parm(string& word, string& valu)
{
  if (!word.compare("HALO_MODEL"))
    switch (atoi(valu.c_str())) {
    case file:
      HALO_MODEL = file;
      break;
    case sing_isothermal:
      HALO_MODEL = sing_isothermal;
      break;
    case isothermal:
      HALO_MODEL = isothermal;
      break;
    case hernquist_model:
      HALO_MODEL = hernquist_model;
      break;
    default:
      cerr << "No such HALO model type: " << (int)HALO_MODEL << endl;
      exit(-2);
    }

  else if (!word.compare("E"))
    E = atof(valu.c_str());

  else if (!word.compare("K"))
    K = atof(valu.c_str());

  else if (!word.compare("INCLINE"))
    INCLINE = atof(valu.c_str());

  else if (!word.compare("PSI"))
    PSI = atof(valu.c_str());

  else if (!word.compare("PHIP"))
    PHIP = atof(valu.c_str());

  else if (!word.compare("VROT"))
    VROT = atof(valu.c_str());

  else if (!word.compare("RCORE"))
    RCORE = atof(valu.c_str());

  else if (!word.compare("RMODMIN"))
    RMODMIN = atof(valu.c_str());

  else if (!word.compare("RMODMAX"))
    RMODMAX = atof(valu.c_str());

  else if (!word.compare("RA"))
    RA = atof(valu.c_str());

  else if (!word.compare("DIVERGE"))
    DIVERGE = atoi(valu.c_str());

  else if (!word.compare("NUMDF"))
    NUMDF = atoi(valu.c_str());

  else if (!word.compare("NRECS"))
    NRECS = atoi(valu.c_str());

  else if (!word.compare("DIVERGE_RFAC"))
    DIVERGE_RFAC = atof(valu.c_str());

  else if (!word.compare("MODFILE"))
    MODFILE = valu;

  else {
    cerr << "No such paramter: " << word << " . . . quitting\n";
    exit(-1);
  }

}

#include <StringTok.H>

string trimLeft(const string &value)
 {
   // stripping only space (character #32)
   string::size_type where = value.find_first_not_of(' ');
  
   if (where == string::npos)
     // string has nothing but space
     return string();
   
   if (where == 0)
     // string has no leading space, don't copy its contents
     return value;
  
   return value.substr(where);
 }


string trimRight(const string &value)
 {
   string::size_type where = value.find_last_not_of(' ');
  
   if (where == string::npos)
     // string has nothing but space
     return string();
   
   if (where == (value.length() - 1))
     // string has no trailing space, don't copy its contents
     return value;
   
   return value.substr(0, where + 1);
 }


string trimComment(const string &value)
 {
   string::size_type where = value.find_first_of("#!");
  
   if (where == string::npos)
     // no comments
     return value;
   
   return value.substr(0, where - 1);
 }


int parse_string(string& line, vector<string>& word, vector<string>& valu)
{
    				// Look for trailing comment
  line = trimComment(line);
  line = trimLeft(trimRight(line));

				// Look for leading . or zero size
  if (line.substr(0,1) == "." || line.size() == 0) return 0;

				// Is there a token marker?
  if (line.find("=") == string::npos) return 0;

				// Parse for tokens
  StringTok<string> tokens(line);
  string first = trimLeft(trimRight(tokens("=")));
  string second = trimLeft(trimRight(tokens("=")));
    
  if (!first.empty()) {
    word.push_back(first);
    valu.push_back(second);
    return 1;
  } 

  return 0;
} 



int get_key_value_from_file(string file, 
			    vector<string>& word, vector<string>& valu)
{
  ifstream in(file.c_str());
  if (!in) {
    cerr << "get_key_value_from_file: error reading <" << file << ">\n";
    return 0;
  }
				// Start with empty elements
  word.erase(word.begin(), word.end());
  valu.erase(word.begin(), word.end());

  int num=0;
  const int lbufsize=1024;
  char lbuf[lbufsize];

  while (in) {
    in.getline(lbuf, lbufsize);
    if (!in) break;
    string line(lbuf);
    num += parse_string(line, word, valu);
  }

  return num;
}


int get_key_value(int argc, char **argv, 
		  vector<string>& word, vector<string>& valu)
{
				// Start with empty elements
  word.erase(word.begin(), word.end());
  valu.erase(word.begin(), word.end());

  int num=0;
  for (int i=0; i<argc; i++) {
    string line(argv[i]);
    num += parse_string(line, word, valu);
  }

  return num;
}




void SatelliteOrbit::parse_args(void)
{
  int i, iret;
  vector<string> word, valu;

  iret = get_key_value_from_file(paramFile, word, valu);
  if (iret == 0) {
    cerr << "SatelliteOrbit: parameter parsing failed\n";
    exit(-1);
  }
				// Set parameters 
  for (i=0; i<iret; i++)
    set_parm(word[i],valu[i]);
}

