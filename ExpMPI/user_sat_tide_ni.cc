
/*
  void user_perturbation()

  Compute acceleration and potential for some desired
  disturbance.

  Turned on by setting flag "user".

*/

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include <math.h>
#include "expand.h"
#undef RNUM			// Conflict

#include <SatelliteOrbit.h>

extern double erf(double);

//  U1	Line of node shift
//  U2	Inclination
//  U3	Time offset for orbit
//  U4	Mass ratio (halo/satellite)
//  U5	Time offset for orbit
//  U6	Turn on duration

SatelliteOrbit* satorb;
static int firstime=1;

extern "C" void user_perturbation()
{
  if (firstime) {
    satorb = new SatelliteOrbit();
				// Convert to radians
    double onedeg = M_PI/180.0;
    double phi   = onedeg*U1;
    double theta = onedeg*U2;
				// Set orientation for satellite disk
    satorb->setTidalOrientation(phi, theta, 0.0);

    firstime = 0;
  }

  double galmass;
  Vector force;

  satorb->setTidalPosition(tnow - U3, 1);

  galmass = U4 * 0.5*(1.0 + erf( (tnow - U5)/U6 ));

  for (int i=1; i<=nbodies; i++)
    {
      if (freeze_particle(i)) continue;

      force = galmass * satorb->tidalForce(x[i] - com2[0], 
					   y[i] - com2[1],
					   z[i] - com2[2],
					   vx[i], vy[i], vz[i]);
		     
      ax[i] += force[1];
      ay[i] += force[2];
      az[i] += force[3];
    }
}

