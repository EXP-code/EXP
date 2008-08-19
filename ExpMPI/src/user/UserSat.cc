#include <math.h>
#include "expand.h"
#include <localmpi.h>
#include <SatelliteOrbit.h>

#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <FindOrb.H>

class UserSat : public ExternalForce
{
private:
  
  string com_name, config, orbfile;
  Component *c0;

  SatelliteOrbit *orb;

  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double core, mass, ton, toff, delta, toffset;

  bool circ, orbit;
  double omega, phase, r0, tlast;

  void userinfo();

public:

				// For debugging . . .
  static int instances;

  UserSat(string &line);
  ~UserSat();

};


int UserSat::instances = 0;

UserSat::UserSat(string &line) : ExternalForce(line)
{

  core = 0.5;			// Satellite core size
  mass = 0.3;			// Satellite mass
  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit

  orbit = false;		// Print out orbit for debugging
  circ = false;			// Assume circular orbit on x-y plane
  r0 = 1.0;			// Radius
  phase = 0.0;			// Initial position angle
  omega = 1.0;			// Angular frequency

  com_name = "sphereSL";	// Default component for com
  config   = "conf.file";	// Configuration file for spherical orbit

  initialize();

				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !com_name.compare(c->name) ) {
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

  if (circ) orb = 0;
  else      orb = new SatelliteOrbit(config);

  if (orbit && myid==0) {
    ostringstream sout;
    sout << "UserSat." << runtag << "." << ++instances;
    orbfile = sout.str();
    ofstream out (orbfile.c_str());
    out << left << setfill('-')
	<< setw(15) << "#"
	<< setw(15) << "+"
	<< setw(15) << "+"
	<< setw(15) << "+"
	<< endl << setfill(' ')
	<< setw(15) << "# Time"
	<< setw(15) << "+ X-pos"
	<< setw(15) << "+ Y-pos"
	<< setw(15) << "+ Z-pos"
	<< endl << setfill('-')
	<< setw(15) << "#"
	<< setw(15) << "+"
	<< setw(15) << "+"
	<< setw(15) << "+"
	<< endl << setfill(' ');
      
    tlast = tnow;
  }

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
  cout << "** User routine SATELLITE IN FIXED POTENTIAL initialized using ";
  if (circ)
    cout << "fixed circular orbit with mass=" << mass 
	 << ", core=" << core
	 << ", r=" << r0 
	 << ", p(0)=" << phase 
	 << ", Omega=" << omega 
	 << endl;
  else
    cout << "specified orbit with mass=" << mass 
	 << ", core=" << core
	 << ", config=" << config
	 << endl;

  print_divider();
}

void UserSat::initialize()
{
  string val;

  if (get_value("comname", val))  com_name = val;
  if (get_value("config", val))   config = val;
  if (get_value("core", val))     core = atof(val.c_str());
  if (get_value("mass", val))     mass = atof(val.c_str());
  if (get_value("ton", val))      ton = atof(val.c_str());
  if (get_value("toff", val))     toff = atof(val.c_str());
  if (get_value("delta", val))    delta = atof(val.c_str());
  if (get_value("toffset", val))  toffset = atof(val.c_str());
  if (get_value("orbit", val))    orbit = atoi(val.c_str()) ? true : false;
  if (get_value("circ", val))     circ = atoi(val.c_str()) ? true : false;
  if (get_value("r0", val))       r0 = atof(val.c_str());
  if (get_value("phase", val))    phase = atof(val.c_str());
  if (get_value("omega", val))    omega = atof(val.c_str());
}

void * UserSat::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], rs[3], fac, ffac, phi;
  double satmass;
  
  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  if (circ) {
    phi = phase + omega*tnow;
    rs[0] = r0*cos(phi);
    rs[1] = r0*sin(phi);
    rs[2] = 0.0;
  }
  else
    orb->get_satellite_orbit(tnow - toffset, &rs[0]);

  satmass = mass * 
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  if (orbit && myid==0 && id==0 && tnow>tlast) {
    ofstream out (orbfile.c_str(), ios::app);
    out << setw(15) << tnow;
    for (int k=0; k<3; k++) out << setw(15) << rs[k];
    out << endl;
    tlast = tnow;
  }


  map<unsigned long, Particle>::iterator it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    fac = core*core;
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k, Component::Inertial) - c0->com[k];
      fac += (pos[k] - rs[k])*(pos[k] - rs[k]);
    }
    fac = pow(fac, -0.5);
    
    ffac = -satmass*fac*fac*fac;

    // Add acceration
    for (int k=0; k<3; k++) cC->AddAcc(i, k, ffac*(pos[k]-rs[k]) );
    
    // Add external potential
    cC->AddPotExt(i, -satmass*fac );
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
