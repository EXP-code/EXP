#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserPeriodic.H>


UserPeriodic::UserPeriodic(string &line) : ExternalForce(line)
{

  id = "PeriodicBC";		// Periodic boundary condition ID

				// Sizes in each dimension
  L = vector<double>(3, 2.0);	
				// Center offset in each dimension
  offset = vector<double>(3, 1.0);

  bc = "ppp";			// Periodic BC in all dimensions;

  comp_name = "";		// Default component

  initialize();
				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !comp_name.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    cerr << "Process " << myid << ": can't find desired component <"
	 << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }
  
  userinfo();
}

UserPeriodic::~UserPeriodic()
{
}

void UserPeriodic::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine PERIODIC BOUNDARY CONDITION initialized, " ;
  cout << ", using component <" << comp_name << ">";
  cout << endl;

  cout << "Cube sides (x , y , x) = (" 
       << L[0] << " , " 
       << L[1] << " , " 
       << L[2] << " ) " << endl; 

  cout << "Center offset (x , y , x) = (" 
       << offset[0] << " , " 
       << offset[1] << " , " 
       << offset[2] << " ) " << endl; 

  cout << "Boundary type (x , y , x) = (" 
       << bc[0] << " , " 
       << bc[1] << " , " 
       << bc[2] << " ) " << endl;

  print_divider();
}

void UserPeriodic::initialize()
{
  string val;

  if (get_value("compname", val))	comp_name = val;

  if (get_value("sx", val))	        L[0] = atof(val.c_str());
  if (get_value("sy", val))	        L[1] = atof(val.c_str());
  if (get_value("sz", val))	        L[2] = atof(val.c_str());

  if (get_value("cx", val))	        offset[0] = atof(val.c_str());
  if (get_value("cy", val))	        offset[1] = atof(val.c_str());
  if (get_value("cz", val))	        offset[2] = atof(val.c_str());

  if (get_value("btype", val)) {
    if (strlen(val.c_str()) >= 3) {
      for (int k=0; k<3; k++) {
	if (val.c_str()[k] == 'p') bc[k] = 'p';
	else                       bc[k] = 'r';
      }
    }
  }
}


void UserPeriodic::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}


void * UserPeriodic::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  
  double pos, delta;
  map<unsigned long, Particle>::iterator it = cC->Particles().begin();
  
  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    
				// Index for the current particle
    unsigned long i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level
    if (multistep && (cC->Part(i)->level < mlevel)) continue;
    
    Particle *p = cC->Part(i);

    for (int k=0; k<3; k++) {

				// Increment so that the positions range
				// between 0 and L[k]
      pos = p->pos[k] + offset[k];

				// Reflection
      if (bc[k] == 'r') {
	if (pos < 0.0) {
	  delta = -pos - L[k]*floor(-pos/L[k]);
	  p->pos[k] = delta - offset[k];
	  p->vel[k] = -p->vel[k];
	} 
	if (pos >= L[k]) {
	  delta = pos - L[k]*floor(pos/L[k]);
	  p->pos[k] =  L[k] - delta - offset[k];
	  p->vel[k] = -p->vel[k];
	}
      }
      
				// Periodic
      if (bc[k] == 'p') {
	if (pos < 0.0) {
	  p->pos[k] = p->pos[k] + L[k]*floor(1.0+fabs(pos/L[k]));
	}
	if (pos >= L[k]) {
	  p->pos[k] = p->pos[k] - L[k]*floor(fabs(pos/L[k]));
	}
      }

    }

    //
    // Sanity check for this particle
    //
    for (int k=0; k<3; k++) {
      if (p->pos[k] < -offset[k] || p->pos[k] >= L[k]-offset[k]) {
	cout << "Process " << myid << " id=" << id 
	     << ": Error in pos[" << k << "]=" << p->pos[k] << endl;
      }
    }
  }
  
  return (NULL);
}


extern "C" {
  ExternalForce *makerPeriodic(string& line)
  {
    return new UserPeriodic(line);
  }
}

class proxypbc { 
public:
  proxypbc()
  {
    factory["userperiodic"] = makerPeriodic;
  }
};

proxypbc p;
