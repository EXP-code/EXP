using namespace std;

#include "expand.h"

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <Vector.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <numerical.h>
#include <orbit.h>
#include <massmodel.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>

#include <HaloBulge.H>

HaloBulge::HaloBulge(string& line) : ExternalForce(line)
{
				// Defaults
  HMODEL = file;
  INFILE = "w05";
  
  MHALO=1.0;
  RHALO=1.0;
  RMODMIN=1.0e-3;
  RMOD=20.0;
  
  RBCORE=1.0;
  MBULGE=1.0;
  RBULGE=1.0;
  RBMODMIN=1.0e-3;
  RBMOD=20.0;

  initialize();

  //
  // Initialize model and orbit
  //

  switch (HMODEL) {
  case file:
    model = new SphericalModelTable(INFILE); // Halo model
    break;
  case isothermal:
    model = new IsothermalSphere(1.0, RMODMIN, RMOD); // Halo model
    break;
  case hernquist_model:
    model = new HernquistSphere(1.0, RMODMIN, RMOD); // Halo model
    break; 
  default:
    cerr << "No such HALO model type: " << (int)HMODEL << endl;
    exit(-2);
  }

  bmodel = new HernquistSphere(RBCORE, RBMODMIN, RBMOD);
}


void * HaloBulge::determine_acceleration_and_potential_thread(void * arg)
{
  double r, potl, dpot, potlB, dpotB;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;
	
    r = 0.0;
    for (int k=0; k<3; k++) r += (*particles)[i].pos[k]*(*particles)[i].pos[k];
    r = sqrt(r);

    model->get_pot_dpot(r/RHALO, potl, dpot);
    potl *= MHALO/RHALO;
    dpot *= MHALO/RHALO/RHALO;

    bmodel->get_pot_dpot(r/RBULGE, potlB, dpotB);
    potlB *= MBULGE/RBULGE;
    dpotB *= MBULGE/RBULGE/RBULGE;

    for (int k=0; k<3; k++) 
      (*particles)[i].acc[k] += -(dpot + dpotB)*(*particles)[i].pos[k]/r;
    (*particles)[i].potext += potl + potlB;
    
  }
}


void HaloBulge::initialize()
{
  string value;

  if (get_value("HMODEL", value))	HMODEL = atoi(value.c_str());
  if (get_value("INFILE", value))	INFILE = value;
  if (get_value("MHALO", value))	MHALO = atof(value.c_str());
  if (get_value("RHALO", value))	RHALO = atof(value.c_str());
  if (get_value("RMODMIN", value))	RMODMIN = atof(value.c_str());
  if (get_value("RMOD", value))		RMOD = atof(value.c_str());
  if (get_value("RBCORE", value))	RBCORE = atof(value.c_str());
  if (get_value("MBULGE", value))	MBULGE = atof(value.c_str());
  if (get_value("RBULGE", value))	RBULGE = atof(value.c_str());
  if (get_value("RBMODMIN", value))	RBMODMIN = atof(value.c_str());
  if (get_value("RBMOD", value))	RBMOD = atof(value.c_str());
}
