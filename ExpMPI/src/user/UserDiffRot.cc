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
  avoid = "";
  maxpm = 2;
  first = true;

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

  if (avoid.size()>0) {
				// Look for component for particle
				// avoidance
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !avoid.compare(c->name) ) {
	c1 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << avoid << ">" << endl;
      c1 = NULL;
    }

  }
  else
    c1 = NULL;

  pos = vector<double>(4*maxpm);
  gen = new ACG(seed);
  uniform = new Uniform(0.0, 1.0, gen);

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
  
  cout << ", avoid = " << avoid;
  cout << ", maxpm = " << maxpm;
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
  if (get_value("avoid", val))		avoid = val;
  if (get_value("maxpm", val))		maxpm = atoi(val.c_str());
  if (get_value("rate", val))		rate = atof(val.c_str());
  if (get_value("seed", val))		seed = atoi(val.c_str());
}


void UserDiffRot::determine_acceleration_and_potential(void)
{
  if (!c0) return;

  // Get particles to avoid
  if (c1) {

    ipm = 0;

    if (myid==0) {

      int number = -1;
      Component::Partstruct *p = c1->get_particles(&number);

      while (p) {
	
	for (int k=0; k<number && ipm<maxpm; k++, ipm++) {
	  pos[ipm*4] = sqrt(p[k].pos[0]*p[k].pos[0] +
			    p[k].pos[1]*p[k].pos[1] +
			    p[k].pos[2]*p[k].pos[2] );
	  pos[ipm*4+1] = p[k].pos[0];
	  pos[ipm*4+2] = p[k].pos[1];
	  pos[ipm*4+3] = p[k].pos[2];
	}

	if (ipm == maxpm) break;
	p = c1->get_particles(&number);
      }
      
      cout << "****************************************************************";
      cout << "UserDiffRot: avoiding " << ipm << " particles\n";
      cout << "****************************************************************";
    }

    MPI_Bcast(&ipm, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&pos[0], ipm*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  exp_thread_fork(false);
}


double UserDiffRot::get_dtime(Particle& p) 
{
  double vv, E, Lx, Ly, Lz, dt;

  // Compute energy
  vv = 0.0;
  for (int k=0; k<3; k++)
    vv += p.vel[k] * p.vel[k];
      
  E = 0.5*vv + p.pot;

  if (E>=0.0) E = -1.0e-10;
  // Compute angular momentum
  Lx = 
    p.pos[1]*p.vel[2]-
    p.pos[2]*p.vel[1];
  
  Ly = 
    p.pos[2]*p.vel[0]-
    p.pos[0]*p.vel[2];
  
  Lz =
    p.pos[0]*p.vel[1]-
    p.pos[1]*p.vel[0];
  
  dt = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)/(-2.0*E);
  
  return dt;
}


void * UserDiffRot::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double phi, cosp, sinp, diffr;
  double xx, yy, uu, vv;

  if (first) {
    
    indx = (*particles)[nbeg].dattrib.size();

    for (int i=nbeg; i<nend; i++)
      (*particles)[i].dattrib.push_back(tnow + 
					get_dtime((*particles)[i])/rate);
    
    first = false;
  }


  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

				// Avoid?
    bool tooclose = false;
    for (int n=0; n<ipm; n++) {
      diffr = 0.0;
      for (int k=0; k<3; k++)
	diffr += 
	  ((*particles)[i].pos[k] - pos[n*4+1+k]) *
	  ((*particles)[i].pos[k] - pos[n*4+1+k]) ;

      if (sqrt(diffr)<pos[n*4]) tooclose = true;
    }
    
    if (tooclose) continue;

    if (tnow>(*particles)[i].dattrib[indx]) {

      (*particles)[i].dattrib[indx] = tnow + get_dtime((*particles)[i])/rate;

				// Do rotation
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
