/*
  First part of leap-frog: step in position
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void incr_position(double dt, int mlevel)
{
  if (eqmotion) {

    list<Component*>::iterator cc;
    Component *c;
  
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;

      unsigned ntot = c->Number();
      for (unsigned i=0; i<ntot; i++) {
				// If we are multistepping, only 
				// advance for this level
      if (multistep && mlevel>=0 && (c->Part(i)->level != mlevel)) continue;

	for (int k=0; k<c->dim; k++) 
	  c->Part(i)->pos[k] += (c->Part(i)->vel[k] - c->covI[k])*dt;
      }
  
      if (c->com_system) {	// Only do this once per multistep
	if (multistep==0 || (mstep==Mstep && mlevel==multistep))
	  for (int k=0; k<c->dim; k++) c->com0[k] += c->cov0[k]*dt;
      }
      
    }

  }

}
