#include <math.h>
#include "expand.h"
#include <SatelliteOrbit.h>

#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPot.H>
#include <biorth.h>
#include <sphereSL.h>
#include <UserResPot.H>
#include <BarForcing.H>


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
  scale = 0.067;
  drfac = 0.05;

  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit
  omega = 18.9;			// Patern speed

  NUME = 400;			// Points in Energy grid
  NUMK = 100;			// Points in Kappa grid

  MASS = 0.05;			// Bar mass
  LENGTH = 0.067;		// Bar length
  COROT = 10;			// Corotation factor
  A21 = 0.2;			// Major to semi-minor ratio
  A32 = 0.05;			// Semi-minor to minor ratio

  domega = 0.0;			// Frequency shift
  t0 = 0.5;			// Mid point of drift
  first = true;

				// Apply force from spherical background
  use_background = true;	// model

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

  initialize();

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;


				// Set up for resonance potential
  SphericalModelTable *hm = new SphericalModelTable(model_file);
  halo_model = hm;

  SphereSL::cache = 0;
  SphereSL::mpi = 1;
  SphereSL *sl = new SphereSL(LMAX, NMAX, NUMR, rmin, rmax, scale, hm);
  halo_ortho = sl;

  ResPot::NUME = NUME;
  ResPot::NUMK = NUMK;
  respot = new ResPot(halo_model, halo_ortho, L0, M0, L1, L2, NMAX);

  BarForcing::L0 = L0;
  BarForcing::M0 = M0;
  BarForcing bar(NMAX, MASS, LENGTH, COROT);
  bar.compute_quad_parameters(A21, A32);
  bar.compute_perturbation(halo_model, halo_ortho, bcoef, bcoefPP);
  omega = bar.Omega();

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
  cout << "** User routine RESONANCE POTENTIAL initialized";
  if (use_background) cout << " using background";
  cout << " with Length=" << LENGTH 
       << ", Mass=" << MASS 
       << ", Omega=" << omega 
       << ", Domega=" << domega 
       << ", T0=" << t0
       << ", b/a=" << A21
       << ", c/b=" << A32
       << ", Ton=" << ton
       << ", Toff=" << toff
       << ", Delta=" << delta
       << "\n";
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
  if (get_value("scale", val))     scale = atof(val.c_str());
  if (get_value("drfac", val))    drfac = atof(val.c_str());

  if (get_value("ton", val))      ton = atof(val.c_str());
  if (get_value("toff", val))     toff = atof(val.c_str());
  if (get_value("delta", val))    delta = atof(val.c_str());
  if (get_value("toffset", val))  toffset = atof(val.c_str());

  if (get_value("MASS", val))     MASS = atof(val.c_str());
  if (get_value("LENGTH", val))   LENGTH = atof(val.c_str());
  if (get_value("COROT", val))    COROT = atof(val.c_str());
  if (get_value("A21", val))      A21 = atof(val.c_str());
  if (get_value("A32", val))      A32 = atof(val.c_str());

  if (get_value("domega", val))   domega = atof(val.c_str());
  if (get_value("t0", val))       t0 = atof(val.c_str());

  if (get_value("NUME", val))     NUME = atoi(val.c_str());
  if (get_value("NUMK", val))     NUMK = atoi(val.c_str());

  if (get_value("use_background", val))
    use_background = atoi(val.c_str()) ? true : false;

  if (get_value("model", val))    model_file = val;
  if (get_value("ctrname", val))  ctr_name = val;
}

void UserResPot::determine_acceleration_and_potential(void)
{
  double Omega = omega*(1.0 + domega*(tnow - t0));
  
  if (first) {
    phase = Omega*tnow;		// Initial phase
    first = false;
  } else {			// Trapezoidal rule integration
    phase += (tnow - tlast)*0.5*(Omega + omlast);
  }

				// Store current state
  tlast = tnow;
  omlast = Omega;

  exp_thread_fork(false);
}

void * UserResPot::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], posN[3], posP[3], vel[3], acc[3];
  double amp, R2, R, pot, dpot;
  
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  amp =
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  for (int i=nbeg; i<nend; i++) {

    R2 = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = (*particles)[i].pos[k];
      if (c0) pos[k] -= c0->com[k];
      vel[k] = (*particles)[i].vel[k];
      R2 += pos[k]*pos[k];
    }
    R = sqrt(R2);

    halo_model->get_pot_dpot(R, pot, dpot);

    for (int k=0; k<3; k++) {

      for (int k2=0; k2<3; k2++) posP[k2] = posN[k2] = pos[k2];

      posN[k] -= drfac*R;
      posP[k] += drfac*R;

      acc[k] = - amp * (
			respot->Pot(posP, vel, phase, bcoef) -
			respot->Pot(posN, vel, phase, bcoef)
			)/(2.0*drfac*R);

      if (use_background && R>0.0) acc[k] += -dpot*pos[k]/R;
    }
    
    for (int k=0; k<3; k++) (*particles)[i].acc[k] += acc[k];
    
    (*particles)[i].potext += amp * 
      respot->Pot(pos, vel, phase, bcoef);

    if (use_background)
      (*particles)[i].potext += (*particles)[i].mass * pot;

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
