/*
  Initialize the velocities of the bodies for the initial timestep
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void init_velocity(void)
{
  int i;

  for (i=1; i<=nbodies; i++) {
    vx[i]=vx[i]+0.5*dtime*ax[i];
    vy[i]=vy[i]+0.5*dtime*ay[i];
    vz[i]=vz[i]+0.5*dtime*az[i];
  }

  /*  Increment velocity time, system time */

  tvel=tvel+0.5*dtime;
  tnow=tvel;
}
