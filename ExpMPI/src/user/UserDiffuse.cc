#include <math.h>
#include <strstream>

#include "expand.h"

#include <ACG.h>
#include <Normal.h>

#include <UserDiffuse.H>

UserDiffuse::UserDiffuse(string &line) : ExternalForce(line)
{
  id = "Energy surface diffusion";

  rate = 0.01;		// Diffusion rate
  name = "";			// Default component name

  initialize();

  if (name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;

  userinfo();
}

UserDiffuse::~UserDiffuse()
{
  delete gen;
  delete nrand;
}

void UserDiffuse::userinfo()
{
  if (myid) return;		// Return if node master node
  if (c0)
    cout << "Diffusion routine initialized for component: <" << name << ">, ";
  else
    cout << "Diffusion routine disabled: no component specified";

  cout << ", rate = " << rate;
  cout << ", seed = " << seed;
  cout << endl;
}

void UserDiffuse::initialize()
{
  string val;

  if (get_value("name", val))		name = val;
  if (get_value("rate", val))		rate = atof(val.c_str());
  if (get_value("seed", val))		seed = atoi(val.c_str());

  gen = new ACG(seed);
  nrand = new Normal(0.0, 1.0, gen);
}


void UserDiffuse::determine_acceleration_and_potential(void)
{
  if (!c0) return;

  exp_thread_fork(false);
}


void * UserDiffuse::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double rr, vv, vv1, dt, dv, sigma, E, vsign;

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

    rr = vv = 0.0;
    for (int k=0; k<3; k++) {
      rr += (*particles)[i].pos[k] * (*particles)[i].pos[k];
      vv += (*particles)[i].vel[k] * (*particles)[i].vel[k];
    }

    rr = sqrt(rr);
    E = 0.5*vv + (*particles)[i].pot;
    if (E>=0.0) continue;

    dt = rr/sqrt(-E);
    sigma = sqrt(fabs(vv*rate)) * dtime/dt;
    
    vv1 = 0.0;
    for (int k=0; k<3; k++) {
     
      dv = sigma*(*nrand)();

      if (dv>0.0) vsign = 1.0;
      else vsign = -1.0;

      dv += (*particles)[i].vel[k];
      
      if (vv1 < vv &&
	  vv1 + dv*dv > vv) dv = vsign*sqrt(vv - vv1);
      else dv = 0.0;
      
      vv1 += dv*dv;

      (*particles)[i].vel[k] = dv;
    }

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerDiffuse(string& line)
  {
    return new UserDiffuse(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["userdiffuse"] = makerDiffuse;
  }
};

proxy p;
