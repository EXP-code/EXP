
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

void satellite_orbit(double T, double* X, double* Y, double* Z);

void user_perturbation()
{
  int i;
   
  double rs[3], fac1, fac2, ffac1, ffac2;
  double satmass;

  satellite_orbit(tnow, &rs[0], &rs[1], &rs[2]);

  satmass = U2 * 0.5*(1.0 + erf( (tnow - U3)/U4 ));

  for (i=1; i<=nbodies; i++)
    {
      if (freeze_particle(i)) continue;

      fac1 = pow(U1*U1 + 
		 (x[i]-rs[0])*(x[i]-rs[0]) +
		 (y[i]-rs[1])*(y[i]-rs[1]) +
		 (z[i]-rs[2])*(z[i]-rs[2]), -0.5);
		     
      fac2 = pow(U1*U1 + 
		 (x[i]+rs[0])*(x[i]+rs[0]) +
		 (y[i]+rs[1])*(y[i]+rs[1]) +
		 (z[i]+rs[2])*(z[i]+rs[2]), -0.5);
		     
      ffac1 = -satmass*fac1*fac1*fac1;
      ffac2 = -satmass*fac2*fac2*fac2;

      ax[i] += ffac1 * (x[i]-rs[0]) + ffac2 * (x[i]+rs[0]);
      ay[i] += ffac1 * (y[i]-rs[1]) + ffac2 * (y[i]+rs[1]);
      az[i] += ffac1 * (z[i]-rs[2]) + ffac2 * (z[i]+rs[2]);
      potext[i] += -satmass*fac1 - satmass*fac2;
    }
}
