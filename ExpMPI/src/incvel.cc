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
  
  list<Component*>::iterator c;
  vector<Particle>::iterator p, pend;
  
  for (c=comp.components.begin(); c != comp.components.end(); c++) {
    
    pend = (*c)->particles.end();
    for (p=(*c)->particles.begin(); p != pend; p++) {
      for (int k=0; k<(*c)->dim; k++) p->vel[k] += p->acc[k]*dtime;
    }
  }

  // Increment times

  tvel=tvel+dtime;
  tnow=tvel;

#ifdef MPE_PROFILE
  MPE_Log_event(14, myid, "e_time");
#endif

}
