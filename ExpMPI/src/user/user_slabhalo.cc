
/*
  void user_perturbation()

  Compute acceleration and potential for some desired
  disturbance.

  Turned on by setting flag "user".

*/

#define HALO

static char rcsid[] = "$Id$";

#include "expand.h"

/*  U1 (0.5)	Satellite core size	*/
/*  U2 (0.3)	Satellite mass		*/
/*  U3 (-20.0)	Turn on start time	*/
/*  U4 (1.0)	Turn on duration	*/
/*  U5 (0.0)	Time offset for orbit   */
/*  U6 (1.0)	Velocity for orbit      */

#ifdef HALO
static double hh = 0.025993723429275473;
static double rho0h = 38.470824928459429;
#endif

static double offset = 0.0;

void user_perturbation()
{
  int i;
#ifdef HALO
  double hfac = 4.0*M_PI*hh*rho0h;
#endif
   
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
#ifdef HALO
      az[i] += ffac * (z[i] - pos2) - hfac*tanh(z[i]/hh);
      potext[i] += -pmass*fac + hfac*hh*log(cosh(z[i]/hh));
#else
      az[i] += ffac * (z[i] - pos2);
      potext[i] += -pmass*fac;
#endif
    }
}
