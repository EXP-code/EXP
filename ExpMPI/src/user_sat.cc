
/*
  void user_perturbation()

  Compute acceleration and potential for some desired
  disturbance.

  Turned on by setting flag "user".

*/

static char rcsid[] = "$Id$";

#include <math.h>
#include "expand.h"
#undef RNUM			// Conflict
#include <SatelliteOrbit.h>

/*  U1 (0.5)	Satellite core size	*/
/*  U2 (0.3)	Satellite mass		*/
/*  U3 (-20.0)	Turn on start time	*/
/*  U4 (1.0)	Turn on duration	*/
/*  U5 (0.0)	Time offset for orbit   */

class UserID {
public:
  UserID() {
    cout << "\n";
    cout << "User routine id:\n";
    cout << "Perturbation using satellite on analytic orbit.\n";
    cout << "\n";
  }
} blab;

SatelliteOrbit *satorb;

extern "C" void user_perturbation()
{
  static bool firstime = true;

  if (firstime) {
    satorb = new SatelliteOrbit;
    firstime = false;
  }

  int i;
   
  double fac, ffac;
  double satmass;

  Vector v = satorb->get_satellite_orbit(tnow - U5);

  satmass = U2 * 0.5*(1.0 + erf( (tnow - U3)/U4 ));

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
