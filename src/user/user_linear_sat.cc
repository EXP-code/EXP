
/*
  void user_perturbation()

  Compute acceleration and potential for some desired
  disturbance.

  Turned on by setting flag "user".

*/

static char rcsid[] = "$Id$";

#include <math.h>
#include <Vector.h>
#include "expand.h"
#undef RNUM			// Conflict

/*  U1 (0.5)	Satellite core size	*/
/*  U2 (0.3)	Satellite mass		*/
/*  U3 (20.0)	Turn on start time	*/
/*  U4 (1.0)	Turn on duration	*/
/*  U5 (1.0)	Satellite velocity      */
/*  U6 (1.0)	Impact parameter        */

extern "C" void user_perturbation()
{
  static bool firstime = true;
  if (firstime && myid==0) {
    cout << "Perturbation using satellite on a linear trajectory\n";
    firstime = false;
  }

  int i;
   
  double fac, ffac;
  double satmass;

  Vector v(1, 3);
  v[1] = U6;
  v[2] = 0.0;
  v[3] = U5*(tnow-U3);

  satmass = U2;

  for (i=1; i<=nbodies; i++)
    {
      if (freeze_particle(i)) continue;

      fac = pow(U1*U1 + 
		(x[i]-v[1])*(x[i]-v[1]) +
		(y[i]-v[2])*(y[i]-v[2]) +
		(z[i]-v[3])*(z[i]-v[3]), -0.5);
		     
      ffac = -satmass*fac*fac*fac;

      ax[i] += ffac * (x[i]-v[1]);
      ay[i] += ffac * (y[i]-v[2]);
      az[i] += ffac * (z[i]-v[3]);
      potext[i] += -satmass*fac;
    }
}
