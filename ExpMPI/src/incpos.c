/*
  First part of leap-frog: step in position
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void incr_position(void)
{
  int i;

#ifdef MPE_PROFILE
  MPE_Log_event(13, myid, "b_time");
#endif

  for(i=1; i<=nbodies; i++) {
    if (freeze_particle(i)) continue;	/* do not move frozen particles
					KL 5/27/92 */
    x[i]=x[i]+vx[i]*dtime;
    y[i]=y[i]+vy[i]*dtime;
    z[i]=z[i]+vz[i]*dtime;
  }

  /* Increment times */

  tpos=tpos+dtime;
  tnow=tpos;

#ifdef MPE_PROFILE
  MPE_Log_event(14, myid, "e_time");
#endif

}
