
/*
  void user_perturbation()

  Compute acceleration and potential for some desired
  disturbance.

  Turned on by setting flag "user".

*/

static char rcsid[] = "$Id$";

#include "expand.h"

#define RS 5.0			/* Satellite galactocentric radius */
#define AS 0.5			/* Satellite core size */
#define MS 0.1			/* Satellite mass */
#define OS 0.2			/* Orbital frequency */
#define TS 10.0			/* Turn on start time */
#define DS 1.0			/* Turn on duration */

void user_perturbation()
{
  int i;
   
  double rs[3], fac, ffac;
  double tfac = 0.5*( 1.0 + erf( (tpos-TS)/DS ));

  rs[0] = RS*cos(tpos*OS);
  rs[1] = RS*cos(tpos*OS);
  rs[2] = 0.0;

  for (i=1; i<=nbodies; i++)
    {
      if (freeze_particle(i)) continue;

      fac = pow(AS*AS + 
		(x[i]-rs[0])*(x[i]-rs[0]) +
		(y[i]-rs[1])*(y[i]-rs[1]) +
		(z[i]-rs[2])*(z[i]-rs[2]), -0.5);
		     
      ffac = -MS*fac*fac*fac;

      ax[i] += ffac * (x[i]-rs[0]);
      ay[i] += ffac * (y[i]-rs[1]);
      az[i] += ffac * (z[i]-rs[2]);
      potext[i] += -MS*fac;
    }
}
