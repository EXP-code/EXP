/*
  Second part of leap-frog: step in velocity
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void incr_velocity(void)
{
#ifdef MPE_PROFILE
  MPE_Log_event(13, myid, "b_time");
#endif
  
  if (eqmotion) {

    list<Component*>::iterator cc;
    vector<Particle>::iterator p, pend;
    unsigned ntot;
    Component *c;
  
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
    
      ntot = c->Number();
      for (int i=0; i<ntot; i++) {

	for (int k=0; k<c->dim; k++) 
	  c->Part(i)->vel[k] += (c->Part(i)->acc[k] - c->acc0[k])*dtime;
      }

      if (c->com_system) {
	for (int k=0; k<c->dim; k++) c->cov0[k] += c->acc0[k]*dtime;
      }
      
    }

  }

  // Increment times

  tvel=tvel+dtime;
  tnow=tvel;

#ifdef MPE_PROFILE
  MPE_Log_event(14, myid, "e_time");
#endif

}
