#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include "expand.h"

#include <tidalField.H>

tidalField::tidalField(string& line) : ExternalForce(line)
{
  hills_omega = 0.5;
  hills_p = 0.5;
  
  initialize();
}

void tidalField::initialize()
{
  string val;

  if (get_value("hills_omega", val)) hills_omega = atof(val.c_str());
  if (get_value("hills_p", val)) hills_p = atof(val.c_str());
}

void * tidalField::determine_acceleration_and_potential_thread(void * arg)
{
  int i;
  double s, c, w2, pp, pm, x, y, z;

  w2 = hills_omega*hills_omega;
  pm = 1.0-hills_p;
  pp = 1.0+hills_p;
	
  c = cos(2.0*hills_omega*tpos);
  s = sin(2.0*hills_omega*tpos);

  int nbeg, nend, id = *((int*)arg);

  list<Component*>::iterator ccp;
  Component *cp;
  for (ccp=comp.components.begin(); ccp != comp.components.end(); ccp++) {
    cp = *ccp;

    nbodies = cp->particles.size();
    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    for (i=nbeg; i<nend; i++) {
      {
	if ((cp->particles)[i].freeze()) continue;
	x = (cp->particles)[i].pos[0];
	y = (cp->particles)[i].pos[1];
	z = (cp->particles)[i].pos[2];
	(cp->particles)[i].acc[0] += 0.5*w2*(pp*(c*x + s*y) - pm*x);
	(cp->particles)[i].acc[1] += 0.5*w2*(pp*(s*x - c*y) - pm*y);
	(cp->particles)[i].acc[2] -= w2*z;
	(cp->particles)[i].potext += 0.5*w2*z*z - 
	  0.25*w2*(pp*(c+s)*x*x + pp*(s-c)*y*y - pm*(x*x+y*y) );
      }
    }
  }

}
