#include <math.h>
#include "expand.h"
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

  double U1, U2, U3, U4, U5;

  void userinfo();

public:

  UserSat(string &line);
  ~UserSat();

};


UserSat::UserSat(string &line) : ExternalForce(line)
{

  U1 = 0.5;			// Satellite core size
  U2 = 0.3;			// Satellite mass
  U3 = -20.0;			// Turn on start time
  U4 = 1.0;			// Turn on duration
  U5  = 0.0;			// Time offset for orbit

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

  if (get_value("comname", val))   com_name = val;
  if (get_value("U1", val))   U1 = atof(val.c_str());
  if (get_value("U2", val))   U2 = atof(val.c_str());
  if (get_value("U3", val))   U3 = atof(val.c_str());
  if (get_value("U4", val))   U4 = atof(val.c_str());
  if (get_value("U5", val))   U5 = atof(val.c_str());
}

void * UserSat::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], rs[3], fac, ffac;
  double satmass;
  
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  orb->get_satellite_orbit(tnow - U5, &rs[0]);

  satmass = U2 * 0.5*(1.0 + erf( (tnow - U3)/U4 ));


  for (int i=nbeg; i<nend; i++) {

    fac = U1*U1;
    for (int k=0; k<3; k++) {
      pos[k] = (*particles)[i].pos[k] - c0->com[k];
      U1 += (pos[k] - rs[k])*(pos[k] - rs[k]);
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

class proxysat { 
public:
  proxysat()
  {
    factory["usersat"] = makerSat;
  }
};

proxysat p;
