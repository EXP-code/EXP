/*
  Synchronize particle coordinates when outputing
  particle data to body data file or when computing diagnostics
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void synchronize_velocity(int sign)
{
  int i;
  double rsign;

  if (sign)
    rsign = -1.0;
  else
    rsign = 1.0;

  for(i=1; i<=nbodies; i++) {
    vx[i]=vx[i]+rsign*ax[i]*0.5*dtime;
    vy[i]=vy[i]+rsign*ay[i]*0.5*dtime;
    vz[i]=vz[i]+rsign*az[i]*0.5*dtime;
  }
  

  /* Increment times */
  tvel=tvel+rsign*0.5*dtime;
  tnow=tvel;

}
