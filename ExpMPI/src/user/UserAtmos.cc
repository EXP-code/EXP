#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserAtmos.H>


UserAtmos::UserAtmos(string &line) : ExternalForce(line)
{

  id = "SphericalHalo";		// Halo model file

				// Components of the acceleration
  g.push_back(-1.0);
  g.push_back( 0.0);
  g.push_back( 0.0);

  compname = "";		// Default component

  initialize();

				// Look for the fiducial component for
				// centering
  c0 = 0;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !compname.compare(c->name) ) {
      c0 = c;
      break;
    }
  }
  
  if (!c0) {
    cerr << "Process " << myid << ": can't find desired component <"
	 << compname << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }


  userinfo();
}

UserAtmos::~UserAtmos()
{
}

void UserAtmos::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine UNIFORM GRAVITATIONAL FIELD initialized" ;
  cout << ", applied to component <" << compname << ">";
  cout << endl;

  cout << "Acceleration constants (gx , gy , gz) = (" 
       << g[0] << " , " 
       << g[1] << " , " 
       << g[2] << " ) " << endl; 

  print_divider();
}

void UserAtmos::initialize()
{
  string val;

  if (get_value("gx", val))	        g[0] = atof(val.c_str());
  if (get_value("gy", val))	        g[1] = atof(val.c_str());
  if (get_value("gz", val))	        g[2] = atof(val.c_str());
  if (get_value("compname", val))	compname = val;
}


void UserAtmos::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}


void * UserAtmos::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  double pot, pos[3];
  
  map<unsigned long, Particle>::iterator it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;


    // Compute potential (up to a some constant)
    pot = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k);	// Inertial by default
      if (c0) pos[k] -= c0->center[k];
      pot += -g[k]*pos[k];
    }

    // Add external accerlation
    for (int k=0; k<3; k++) cC->AddAcc(i, k, g[k]);

    // Add external potential
    cC->AddPotExt(i, pot);

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerAtmos(string& line)
  {
    return new UserAtmos(line);
  }
}

class proxyatmos { 
public:
  proxyatmos()
  {
    factory["useratmos"] = makerAtmos;
  }
};

proxyatmos p;
