#include <math.h>
#include <strstream>

#include "expand.h"

#include <ACG.h>
#include <Normal.h>

#include <UserDiffRot.H>

UserDiffRot::UserDiffRot(string &line) : ExternalForce(line)
{
  id = "Rotational randomization";

  rate = 0.5;			// Rate relative to dyn time
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

UserDiffRot::~UserDiffRot()
{
  delete gen;
  delete uniform;
}

void UserDiffRot::userinfo()
{
  if (myid) return;		// Return if node master node
  if (c0)
    cout << "** User routine ROTATION RANDOMIZATION initialized for component: <" << name << ">";
  else
    cout << "** User routine ROTATION RANDOMIZATION disabled: no component specified";
  
  cout << ", rate = " << rate;
  cout << ", seed = " << seed;
  cout << endl;
  cout << "****************************************************************"
       << endl;

}

void UserDiffRot::initialize()
{
  string val;

  if (get_value("name", val))		name = val;
  if (get_value("rate", val))		rate = atof(val.c_str());
  if (get_value("seed", val))		seed = atoi(val.c_str());

  gen = new ACG(seed);
  uniform = new Uniform(0.0, 1.0, gen);
}


void UserDiffRot::determine_acceleration_and_potential(void)
{
  if (!c0) return;

  exp_thread_fork(false);
}


void * UserDiffRot::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double E, dt, phi, cosp, sinp;
  double Lx, Ly, Lz, xx, yy, uu, vv;

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

				// Compute energy
    vv = 0.0;
    for (int k=0; k<3; k++)
      vv += (*particles)[i].vel[k] * (*particles)[i].vel[k];

    E = 0.5*vv + (*particles)[i].pot;

    if (E>=0.0) continue;

				// Compute angular momentum
    Lx = 
      (*particles)[i].pos[1]*(*particles)[i].vel[2]-
      (*particles)[i].pos[2]*(*particles)[i].vel[1];

    Ly = 
      (*particles)[i].pos[2]*(*particles)[i].vel[0]-
      (*particles)[i].pos[0]*(*particles)[i].vel[2];

    Lz =
      (*particles)[i].pos[0]*(*particles)[i].vel[1]-
      (*particles)[i].pos[1]*(*particles)[i].vel[0];

    dt = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)/sqrt(-2.0*E);
    
				// Do the rotation
    if ((*uniform)() < dtime*rate/dt) {
      phi = 2.0*M_PI*(*uniform)();
      cosp = cos(phi);
      sinp = sin(phi);
      
      xx = (*particles)[i].pos[0]*cosp - (*particles)[i].pos[1]*sinp;
      yy = (*particles)[i].pos[0]*sinp + (*particles)[i].pos[1]*cosp;
      uu = (*particles)[i].vel[0]*cosp - (*particles)[i].vel[1]*sinp;
      vv = (*particles)[i].vel[0]*sinp + (*particles)[i].vel[1]*cosp;

      (*particles)[i].pos[0] = xx;
      (*particles)[i].pos[1] = yy;
      (*particles)[i].vel[0] = uu;
      (*particles)[i].vel[1] = vv;
    }

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerDiffRot(string& line)
  {
    return new UserDiffRot(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["userdiffrot"] = makerDiffRot;
  }
};

proxy p;
