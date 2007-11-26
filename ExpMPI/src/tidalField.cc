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
	
  c = cos(2.0*hills_omega*tnow);
  s = sin(2.0*hills_omega*tnow);

  int nbeg, nend, id = *((int*)arg);

  list<Component*>::iterator ccp;
  Component *cp;
  for (ccp=comp.components.begin(); ccp != comp.components.end(); ccp++) {
    cp = *ccp;

    nbodies = cp->Number();
    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    map<unsigned long, Particle>::iterator it = cp->Particles().begin();
    unsigned long i;

    for (int q=0   ; q<nbeg; q++) it++;
    for (int q=nbeg; q<nend; q++) 
      {
	i = (it++)->first;

	if (cp->freeze(i)) continue;
	x = cp->Pos(i, 0);
	y = cp->Pos(i, 1);
	z = cp->Pos(i, 2);
	cp->AddAcc(i, 0, 0.5*w2*(pp*(c*x + s*y) - pm*x) );
	cp->AddAcc(i, 1, 0.5*w2*(pp*(s*x - c*y) - pm*y) );
	cp->AddAcc(i, 2, w2*z );
	cp->AddPotExt(i, 0.5*w2*z*z - 
		      0.25*w2*(pp*(c+s)*x*x + pp*(s-c)*y*y - pm*(x*x+y*y) ) );
      }
  }
}
