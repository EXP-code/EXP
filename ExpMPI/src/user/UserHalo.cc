#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserHalo.H>

UserHalo::UserHalo(string &line) : ExternalForce(line)
{

  id = "SphericalHalo";		// Halo model file
  model_file = "SLGridSph.model";

  diverge = 0;			// Use analytic divergence (true/false)
  diverge_rfac = 1.0;		// Exponent for profile divergence
  ctr_name = "";		// Default component for center

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


  model = new SphericalModelTable(model_file, diverge, diverge_rfac);

  userinfo();
}

UserHalo::~UserHalo()
{
  delete model;
}

void UserHalo::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SPHERICAL HALO initialized, " ;
  cout << "Filename=" << model_file << "  diverge=" << diverge
       << "  diverge_rfac=" << diverge_rfac;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << endl;
  
  print_divider();
}

void UserHalo::initialize()
{
  string val;

  if (get_value("model_file", val))	model_file = val;
  if (get_value("diverge", val))	diverge = atoi(val.c_str());
  if (get_value("diverge_rfac", val))	diverge_rfac = atof(val.c_str());
  if (get_value("ctrname", val))	ctr_name = val;
}


void UserHalo::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}


void * UserHalo::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double xx, yy, zz, rr, r, pot, dpot, pos[3];

  for (int i=nbeg; i<nend; i++) {

    rr = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = (*particles)[i].pos[k];
      if (c0) pos[k] -= c0->center[k];
      rr += pos[k]*pos[i];
    }
    r = sqrt(rr);

    model->get_pot_dpot(r, pot, dpot);

    for (int k=0; k<3; k++)
      (*particles)[i].acc[k] += -dpot*pos[k]/r;
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerHalo(string& line)
  {
    return new UserHalo(line);
  }
}

class proxyhalo { 
public:
  proxyhalo()
  {
    factory["userhalo"] = makerHalo;
  }
};

proxyhalo p;
