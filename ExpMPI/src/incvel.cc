/*
  Second part of leap-frog: step in velocity
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void incr_velocity(double dt, int mlevel)
{
  if (eqmotion) {

    list<Component*>::iterator cc;
    vector<Particle>::iterator p, pend;
    Component *c;
  
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
    
      unsigned ntot = c->Number();
      for (unsigned i=0; i<ntot; i++) {
				// If we are multistepping, only 
				// advance for this level
	if (multistep && (c->Part(i)->level != mlevel)) continue;

	for (int k=0; k<c->dim; k++) 
	  c->Part(i)->vel[k] += (c->Part(i)->acc[k] - c->acc0[k])*dt;
      }

      if (c->com_system) {	// Only do this once per multistep
	if (multistep==0 || (mstep==Mstep && mlevel==multistep))
	  for (int k=0; k<c->dim; k++) c->cov0[k] += c->acc0[k]*dt;
      }
      
    }

  }

}
