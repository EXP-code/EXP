
/*
  void tidal_field()

  	Compute tidal acceleration and potential in the Hills limit.
  The calculation is done in a frame having its cartesian axis 
  directions fixed in space; i.e., the system is *not* phase
  locked with the perturber as is usually the case.
  Unfortunately, I know of no conserved quantities in this frame.


  int freeze_particle(int i)

	Freeze the i-th particle if it has gone too far away from the
	origin.


  KL 5/27/92

*/

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#include "expand.h"

void tidal_field()
{
	int i;
	double s, c, w2, pp, pm;

	w2 = hills_omega*hills_omega;
	pm = 1.0-hills_p;
	pp = 1.0+hills_p;

	c = cos(2.0*hills_omega*tpos);
	s = sin(2.0*hills_omega*tpos);

	for (i=1; i<=nbodies; i++)
	{
		if (freeze_particle(i)) continue;
		ax[i] += 0.5*w2*(pp*(c*x[i] + s*y[i]) - pm*x[i]);
		ay[i] += 0.5*w2*(pp*(s*x[i] - c*y[i]) - pm*y[i]);
		az[i] -= w2*z[i];
		potext[i] += 0.5*w2*z[i]*z[i] - 0.25*w2*(pp*(c+s)*x[i]*x[i]
			+ pp*(s-c)*y[i]*y[i] - pm*(x[i]*x[i]+y[i]*y[i]) );
	}
}




/*

int freeze_particle(int i)
{
	double v2, pe, ke;

	v2 = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	ke = 0.5*v2;
	pe = pot[i];
	if (pe != 0.0 && ke + pe > 0.0) return 1;
	else return 0;
}



*/
	


int freeze_particle(int i)
{
	double r2;


	r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
	if (r2 > rmax_tidal*rmax_tidal) 
	{
		x[i] = 1.e10;
		y[i] = 1.e10;
		z[i] = 1.e10;
		return 1;
	}
	else return 0;
}



	
				/* Compute shock perturbation */
				/* Added by MDW 3/25/94 */

void external_shock()
{
  int i;
  double w2;
  double get_tidal_shock(double);


  w2 = get_tidal_shock(tpos);

  for (i=1; i<=nbodies; i++)
    {
      if (freeze_particle(i)) continue;
      az[i] -= w2*z[i];
      potext[i] += 0.5*w2*z[i]*z[i];
    }
}
