/*
  Second part of leap-frog: step in velocity
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void incr_velocity(void)
{
  int i;

#ifdef MPE_PROFILE
  MPE_Log_event(13, myid, "b_time");
#endif

  if (myid != 0) {

    for (i=1; i<=nbodies; i++) {
      vx[i]=vx[i]+ax[i]*dtime;
      vy[i]=vy[i]+ay[i]*dtime;
      vz[i]=vz[i]+az[i]*dtime;
    }
  }

  /* Increment times */

  tvel=tvel+dtime;
  tnow=tvel;

#ifdef MPE_PROFILE
  MPE_Log_event(14, myid, "e_time");
#endif

}
