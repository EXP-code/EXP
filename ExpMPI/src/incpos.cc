/*
  First part of leap-frog: step in position
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void incr_position(void)
{
#ifdef MPE_PROFILE
  MPE_Log_event(13, myid, "b_time");
#endif

  if (eqmotion) {

    list<Component*>::iterator cc;
    vector<Particle>::iterator p, pend;
    Component *c;
  
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;

      pend = c->particles.end();
      for (p=c->particles.begin(); p != pend; p++) {
	
	if (c->freeze(*p)) continue;
	
	for (int k=0; k<c->dim; k++) 
	  p->pos[k] += (p->vel[k] - c->cov0[k])*dtime;
      }
  
      if (c->com_system) {
	for (int k=0; k<c->dim; k++) c->com0[k] += c->cov0[k]*dtime;
      }
      
    }

  }

  // Increment times

  tpos=tpos+dtime;
  tnow=tpos;

#ifdef MPE_PROFILE
  MPE_Log_event(14, myid, "e_time");
#endif

}
