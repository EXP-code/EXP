#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserSlabHalo.H>


UserSlabHalo::UserSlabHalo(string &line) : ExternalForce(line)
{

  id   = "SlabHalo";		// Halo model file

  rho0 = 1.0;			// Central density
  h0   = 1.0;			// Scale height
  U0   = 4.0*M_PI*rho0*h0*h0;	// Potential coefficient

				// Midplane squared velocity for 
  v20  = 2.0*U0*log(cosh(1.0));	// turning point at h0

  v0   = sqrt(v20);	        // Midplane velocity for turning point at h0

  ctr_name = "";		// Default component for center (none)

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
      cerr << "libslabhalo, process " << myid 
	   << ": can't find desired component <"
	   << ctr_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;


  userinfo();
}

UserSlabHalo::~UserSlabHalo()
{
}

void UserSlabHalo::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SLAB HALO initialized, v0=" << v0 
       << ", rho0=" << rho0 << ", h0=" << h0;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << endl;

  print_divider();
}

void UserSlabHalo::initialize()
{
  string val;

  if (get_value("ctrname", val))	ctr_name = val;

  if (get_value("h0", val))	        h0 = atof(val.c_str());

  if (get_value("rho0", val)) {
    rho0 = atof(val.c_str());
    U0  = 4.0*M_PI*rho0*h0*h0;
    v20 = 2.0*U0*log(cosh(1.0));
    v0 = sqrt(v20);
  }

  if (get_value("v0", val)) {
    v0 = atof(val.c_str());
    v20 = v0*v0;
    U0 = 0.5*v20/log(cosh(1.0));
    rho0 = U0/(4.0*M_PI*h0*h0);
  }

}


void UserSlabHalo::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);

  print_timings("UserSlabHalo: accleration timings");
}


void * UserSlabHalo::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  double pos[3];

  map<unsigned long, Particle>::iterator it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level
    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k);	// Inertial by default
      if (c0) pos[k] -= c0->center[k];
    }
    
    // Add external accerlation
    cC->AddAcc(i, 2, -U0/h0*tanh(pos[2]/h0));
    
    // Add external potential
    cC->AddPotExt(i, U0*log(cosh(pos[2]/h0)) );
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerSlabHalo(string& line)
  {
    return new UserSlabHalo(line);
  }
}

class proxyslabhalo { 
public:
  proxyslabhalo()
  {
    factory["userslabhalo"] = makerSlabHalo;
  }
};

proxyslabhalo p;
