
/*
  void user_perturbation()

  Compute acceleration and potential for some desired
  disturbance.

  Turned on by setting flag "user".

*/

static char rcsid[] = "$Id$";

#include "expand.h"

/*  U1 (0.5)	Satellite core size	*/
/*  U2 (0.3)	Satellite mass		*/
/*  U3 (-20.0)	Turn on start time	*/
/*  U4 (1.0)	Turn on duration	*/
/*  U5 (0.0)	Time offset for orbit   */
/*  U6 (1.0)	Velocity for orbit      */

static double offset = 0.0;

void user_perturbation()
{
  int i;
   
  double rs[3], fac, ffac;
  double pmass, pos2;

  offset -= zcm_slab;

  pos2 = U6*(tnow - U5) + offset;
  
  pmass = U2 * 0.5*(1.0 + erf( (tnow - U3)/U4 ));

  for (i=1; i<=nbodies; i++)
    {
      if (freeze_particle(i)) continue;

      fac = pow(U1*U1 + 
		(x[i] - 0.5)*(x[i] - 0.5) + 
		(y[i] - 0.5)*(y[i] - 0.5) + 
		(z[i] - pos2)*(z[i] - pos2), -0.5);
		     
      ffac = -pmass*fac*fac*fac;

      ax[i] += ffac * (x[i] - 0.5);
      ay[i] += ffac * (y[i] - 0.5);
      az[i] += ffac * (z[i] - pos2);
      potext[i] += -pmass*fac;
    }
}
