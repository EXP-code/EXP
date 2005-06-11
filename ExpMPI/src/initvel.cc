/*
  Initialize the velocities of the bodies for the initial timestep
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void init_velocity(void)
{
  list<Component*>::iterator cc;
  vector<Particle>::iterator p, pend;
  unsigned ntot;
  Component *c;

  if (eqmotion) {

    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;

      ntot = c->Number();
      for (int i=0; i<ntot; i++) {

	for (int k=0; k<c->dim; k++) 
	  c->Part(i)->vel[k] += 0.5*dtime*c->Part(i)->acc[k];
      }
      
      if (c->com_system) {
	for (int k=0; k<c->dim; k++) c->cov0[k] += 0.5*dtime*c->acc0[k];
      }
      
    }
    
  }

  // Increment velocity time, system time

  tvel=tvel+0.5*dtime;
  tnow=tvel;
}
