#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserShear.H>


UserShear::UserShear(string &line) : ExternalForce(line)
{

  id = "ShearingSheet";		// Shearing sheet
  r0 = 1;			// R_o in radial units
  s0 = 1;			// Sigma in velocity units
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

  userinfo();
}

UserShear::~UserShear()
{
}

void UserShear::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SHEARING SHEET initialized, " ;
  cout << "r_0=" << r0 << "  sigma=" << s0 << "  v_0=" << sqrt(2.0)*s0;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << endl;

  print_divider();
}

void UserShear::initialize()
{
  string val;

  if (get_value("radius", val))	        r0 = atof(val.c_str());
  if (get_value("velocity", val))       s0 = atof(val.c_str());
  if (get_value("ctrname", val))	ctr_name = val;

  omega = sqrt(2.0)*s0/r0;
  kappa = 2.0*s0/r0;
  A     =  s0/(sqrt(2.0)*r0);
  B     = -s0/(sqrt(2.0)*r0);
}


void UserShear::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);

  print_timings("UserShear: accleration timings");
}


void * UserShear::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  double x, xd, yd, ax, ay, pot;

  PartMapItr it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    unsigned long i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    x  = cC->Pos(i, 0);
    xd = cC->Vel(i, 0);
    yd = cC->Vel(i, 1);

    ax =  2.0*omega*yd + 4.0*A*omega*x;
    ay = -2.0*omega*xd;

    cC->AddAcc(i, 0, ax);
    cC->AddAcc(i, 0, ay);

    // Add external potential
    pot = -2.0*A*omega*x*x;
    cC->AddPotExt(i, pot);
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerShear(string& line)
  {
    return new UserShear(line);
  }
}

class proxyhalo { 
public:
  proxyhalo()
  {
    factory["usershear"] = makerShear;
  }
};

proxyhalo p;
