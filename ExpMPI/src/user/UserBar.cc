#include <math.h>
#include "expand.h"

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

class UserBar : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double U1, U2, U3, U4, U5;

  void userinfo();

public:

  UserBar(string &line);
  ~UserBar();

};


UserBar::UserBar(string &line) : ExternalForce(line)
{
  id = "RotatingBar";

  U1 = 0.5;			// Bar length
  U2 = 0.3;			// Bar amplitude
  U3 = -20.0;			// Turn on start time
  U4 = 1.0;			// Turn on duration
  U5  = 0.0;			// Corotation factor

  com_name = "";		// Default component for com

  initialize();

  if (com_name.size()>0) {
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

  }
  else
    c0 = NULL;

  userinfo();
}

UserBar::~UserBar()
{
}

void UserBar::userinfo()
{
  if (myid) return;		// Return if node master node
  cout << "User routine initialized: rotating bar with fixed pattern speed, " ;
  if (c0) 
    cout << " center on component <" << com_name << ">" << endl;
  else
    cout << " center on origin" << endl;
}

void UserBar::initialize()
{
  string val;

  if (get_value("comname", val))   com_name = val;
  if (get_value("U1", val))   U1 = atof(val.c_str());
  if (get_value("U2", val))   U2 = atof(val.c_str());
  if (get_value("U3", val))   U3 = atof(val.c_str());
  if (get_value("U4", val))   U4 = atof(val.c_str());
  if (get_value("U5", val))   U5 = atof(val.c_str());
}

void * UserBar::determine_acceleration_and_potential_thread(void * arg) 
{
  static const double numfac = 3.86274202023190e-01;
  static bool firstime = true;
  static double posang, lastomega, lasttime;

  list<Component*>::iterator cc;
  Component *c;

				// Get frequency
  double R=U1*U5, theta=1.57079632679490, phi;
  double dens, potl, potr, pott, potp;
  double avg = 0.0;
  
  for (int n=0; n<8; n++) {
    phi = 2.0*M_PI/8.0 * n;

    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;

      if (c->force->geometry == PotAccel::sphere || 
	  c->force->geometry == PotAccel::cylinder) {

	((Basis*)c->force)->
	  determine_fields_at_point_sph(R, theta, phi,
					&dens, &potl, &potr, &pott, &potp);

	avg += potr/8.0;
      }
    }
  }

  double omega = sqrt(avg/R);

  if (firstime) {
    posang = 0.0;
    lastomega = omega;
    lasttime = tnow;
    firstime = false;
  } else {
    posang += 0.5*(omega + lastomega)*dtime;
    lastomega = omega;
    lasttime = tnow;
#ifdef DEBUG
    if (myid==0)
      cout << "Time=" << tnow << " Posang=" << posang << endl;
#endif
  }


  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double fac, ffac, amp = U2 * 0.5*(1.0 + erf( (tnow - U3)/U4 ));
  double xx, yy, zz, rr, nn;
  vector<double> pos(3); 
  double cos2p = cos(2.0*posang);
  double sin2p = sin(2.0*posang);

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

    if (c0)
      for (int k=0; k<3; k++) pos[k] = (*particles)[i].pos[k] - c0->com[k];
    else
      for (int k=0; k<3; k++) pos[k] = (*particles)[i].pos[k];
    
    xx = pos[0];
    yy = pos[1];
    zz = pos[2];
    rr = sqrt( xx*xx + yy*yy + zz*zz );

    fac = U1 + rr;
		     
    ffac = -amp*numfac*U1*U1*U1/pow(fac, 6.0);

    nn = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
      
    (*particles)[i].acc[0] += ffac*
      (  2.0*xx*cos2p*fac + 2.0*yy*sin2p*fac - 5.0*nn*xx/(rr+1.0e-16) );
    
    (*particles)[i].acc[1] += ffac*
      ( -2.0*yy*cos2p*fac + 2.0*xx*sin2p*fac - 5.0*nn*yy/(rr+1.0e-16) );

    (*particles)[i].acc[2] += ffac*
      ( -5.0*nn*zz/(rr+1.0e-16) );
    
    (*particles)[i].potext += -ffac*nn*fac;
    
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerBar(string& line)
  {
    return new UserBar(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["userbar"] = makerBar;
  }
};

proxy p;
