#include <math.h>
#include "expand.h"
#include <SatelliteOrbit.h>

#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPot.H>
#include <biorth.h>
#include <biorthSL.h>
#include <UserResPot.H>
#include <BarForcing.H>


SphericalModelTable *hm;
                                                                                
double pot(double r)
{
  return hm->get_pot(r);
}
                                                                                
double dpot(double r)
{
  return hm->get_dpot(r);
}
                                                                                
double dens(double r)
{
  return hm->get_density(r);
}

UserResPot::UserResPot(string &line) : ExternalForce(line)
{
  LMAX = 2;
  NMAX = 20;
  NUMR = 1000;
  L0 = 2;
  M0 = 2;
  L1 = -1;
  L2 =  2;
  rmin = 1.0e-4;
  rmax = 1.95;
  drfac = 0.05;

  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit
  omega = 18.9;			// Patern speed

  MASS = 0.05;			// Bar mass
  LENGTH = 0.067;		// Bar length
  COROT = 10;			// Corotation factor

				// Tabled spherical model
  model_file = "SLGridSph.model";
  com_name = "sphereSL";	// Default component for com

  initialize();

				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !com_name.compare(c->id) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    cerr << "Process " << myid << ": can't find desired component <"
	 << com_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }


  hm =  new SphericalModelTable(model_file);
  halo_model = hm;

  SphereSL::cache = 0;
  halo_ortho = new SphereSL(LMAX, NMAX, NUMR, rmin, rmax,
			    pot, dpot, dens);

  respot = new ResPot(halo_model, halo_ortho, L0, M0, L1, L2, NMAX);

  BarForcing::L0 = L0;
  BarForcing::M0 = M0;
  BarForcing bar(NMAX, MASS, LENGTH, COROT);
  bar.compute_quad_parameters();
  bar.compute_perturbation(halo_model, halo_ortho, bcoef, bcoefPP);

  userinfo();
}

UserResPot::~UserResPot()
{
  delete halo_model;
  delete halo_ortho;
}

void UserResPot::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine SATELLITE IN FIXED POTENTIAL initialized\n";
  print_divider();
}

void UserResPot::initialize()
{
  string val;

  if (get_value("LMAX", val))     LMAX = atoi(val.c_str());
  if (get_value("NMAX", val))     NMAX = atoi(val.c_str());
  if (get_value("NUMR", val))     NUMR = atoi(val.c_str());

  if (get_value("L0", val))       L0 = atoi(val.c_str());
  if (get_value("M0", val))       M0 = atoi(val.c_str());
  if (get_value("L1", val))       L1 = atoi(val.c_str());
  if (get_value("L2", val))       L2 = atoi(val.c_str());

  if (get_value("rmin", val))     rmin = atof(val.c_str());
  if (get_value("rmax", val))     rmax = atof(val.c_str());
  if (get_value("drfac", val))    drfac = atof(val.c_str());

  if (get_value("ton", val))      ton = atof(val.c_str());
  if (get_value("toff", val))     toff = atof(val.c_str());
  if (get_value("delta", val))    delta = atof(val.c_str());
  if (get_value("toffset", val))  toffset = atof(val.c_str());
  if (get_value("omega", val))    omega = atof(val.c_str());

  if (get_value("MASS", val))     MASS = atof(val.c_str());
  if (get_value("LENGTH", val))   LENGTH = atof(val.c_str());
  if (get_value("COROT", val))    COROT = atof(val.c_str());


  if (get_value("file", val))     model_file = val;
}

void * UserResPot::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], posN[3], posP[3], vel[3], acc[3];
  double amp, R2, R=0.0;
  
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  amp =
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    

  for (int i=nbeg; i<nend; i++) {

    R2 = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = (*particles)[i].pos[k] - c0->com[k];
      vel[k] = (*particles)[i].vel[k];
      R2 += pos[k]*pos[k];
    }

    for (int k=0; k<3; k++) {

      for (int k2=0; k<3; k++) posP[k2] = posN[k2] = pos[k2];

      posN[k] -= drfac*R;
      posP[k] += drfac*R;

      acc[k] = - amp * (
			respot->Pot(posP, vel, omega, tpos, bcoef) -
			respot->Pot(posN, vel, omega, tpos, bcoef)
			)/(2.0*drfac*R);
    }
    
    for (int k=0; k<3; k++) (*particles)[i].acc[k] += acc[k];
    
    (*particles)[i].potext += amp * 
      respot->Pot(pos, vel, omega, tpos, bcoef);
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerResPot(string& line)
  {
    return new UserResPot(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["userrespot"] = makerResPot;
  }
};

static proxy p;
