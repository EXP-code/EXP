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
    Component *c;
    int ntot;
  
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;

      ntot = c->Number();
      for (int i=0; i<ntot; i++) {
	
	for (int k=0; k<c->dim; k++) 
	  c->Part(i)->pos[k] += (c->Part(i)->vel[k] - c->covI[k])*dtime;
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
