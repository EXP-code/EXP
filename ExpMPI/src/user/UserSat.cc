#include <math.h>
#include "expand.h"
#include <localmpi.h>
#include <SatelliteOrbit.h>

#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

class UserSat : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  SatelliteOrbit *orb;

  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double core, mass, ton, toff, delta, toffset;

  void userinfo();

public:

  UserSat(string &line);
  ~UserSat();

};


UserSat::UserSat(string &line) : ExternalForce(line)
{

  core = 0.5;			// Satellite core size
  mass = 0.3;			// Satellite mass
  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit

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

  orb = new SatelliteOrbit;

  userinfo();
}

UserSat::~UserSat()
{
  delete orb;
}

void UserSat::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine SATELLITE IN FIXED POTENTIAL initialized\n";
  print_divider();
}

void UserSat::initialize()
{
  string val;

  if (get_value("comname", val))  com_name = val;
  if (get_value("core", val))     core = atof(val.c_str());
  if (get_value("mass", val))     mass = atof(val.c_str());
  if (get_value("ton", val))      ton = atof(val.c_str());
  if (get_value("toff", val))     toff = atof(val.c_str());
  if (get_value("delta", val))    delta = atof(val.c_str());
  if (get_value("toffset", val))  toffset = atof(val.c_str());
}

void * UserSat::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], rs[3], fac, ffac;
  double satmass;
  
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  orb->get_satellite_orbit(tnow - toffset, &rs[0]);

  satmass = mass * 
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    


  for (int i=nbeg; i<nend; i++) {

    fac = core*core;
    for (int k=0; k<3; k++) {
      pos[k] = (*particles)[i].pos[k] - c0->com[k];
      fac += (pos[k] - rs[k])*(pos[k] - rs[k]);
    }
    fac = pow(fac, -0.5);
    
    ffac = -satmass*fac*fac*fac;

    for (int k=0; k<3; k++)
      (*particles)[i].acc[k] += ffac*(pos[k]-rs[k]);
    
    (*particles)[i].potext += -satmass*fac;
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerSat(string& line)
  {
    return new UserSat(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["usersat"] = makerSat;
  }
};

static proxy p;
