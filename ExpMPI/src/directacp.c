/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the direct summation
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *      KL 5/27/92   Modified to allow freezing of particles beyond some cutoff
 *    radius. This is needed when running systems in tidal fields.
 *    The function freeze_particle() is called for each particle
 *    to decide whether to freeze.
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *
 ***************************************************************************/

#include "expand.h"

static char rcsid[] = "$Id$";

void determine_acceleration_and_potential_n2(void);

static double soft2;		/* Particle softening squared */

void get_acceleration_and_potential_n2(void)
{
  static int firstime=1;

  if (firstime) {

    /* Initialize softening */

    soft2 = soft*soft;
    firstime = 0;
  }


  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  determine_acceleration_and_potential_n2();

}

void determine_acceleration_and_potential_n2(void)
{
  int i, j;
  double rr, rfac;

  /* Determine potential and acceleration */

  for (i=1; i<=nbodies; i++) {

    ax[i] = 0.0;
    ay[i] = 0.0;
    az[i] = 0.0;
    pot[i] = 0.0;

    if (component[i] != 1) continue;

    if (freeze_particle(i)) 	/* frozen particles do not respond */
    {				/* KL 5/27/92 */
	continue;
    }

    for (j=1; j<i; j++) {

      rr = sqrt( (x[i]-x[j])*(x[i]-x[j]) +
		 (y[i]-y[j])*(y[i]-y[j]) +
		 (z[i]-z[j])*(z[i]-z[j]) + soft2);

      rfac = 1.0/(rr*rr*rr);

      ax[i] += -mass[j]*(x[i] - x[j])*rfac;
      ay[i] += -mass[j]*(y[i] - y[j])*rfac;
      az[i] += -mass[j]*(z[i] - z[j])*rfac;
      pot[i] += -mass[j]/rr;
    }

    for (j=i+1; j<=nbodies; j++) {

      rr = sqrt( (x[i]-x[j])*(x[i]-x[j]) +
		 (y[i]-y[j])*(y[i]-y[j]) +
		 (z[i]-z[j])*(z[i]-z[j]) + soft2);

      rfac = 1.0/(rr*rr*rr);

      ax[i] += -mass[j]*(x[i] - x[j])*rfac;
      ay[i] += -mass[j]*(y[i] - y[j])*rfac;
      az[i] += -mass[j]*(z[i] - z[j])*rfac;
      pot[i] += -mass[j]/rr;
    }

  }
  
}

