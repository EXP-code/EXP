/* 
  Call necessary routines to advance phase-space one step
*/

#include "expand.h"

void adjust_scale(int);
void out_put(int);

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void do_step(int n)
{

  /*========================*/
  /* Advance using leapfrog */
  /*========================*/

				/* put COM at origin in phase space */
  if (fixpos==1) fix_positions();
  if (fixpos==2) fix_positions_by_component();

				/* Position step */
  incr_position();

				/* Compute acceleration */
  compute_potential();

				/* Adjust to keep COM at origin */
  if (fixacc) fix_acceleration();

				/* Velocity step */
  incr_velocity();

				/* Adjust scale length */

  if (adjust) adjust_scale(n);


/*
  if (myid != 0)
  fprintf(stderr, "Process %d: %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", 
	  myid,
	  x[nbodies/2], y[nbodies/2], z[nbodies/2], 
	  vx[nbodies/2], vy[nbodies/2], vz[nbodies/2],
	  ax[nbodies/2], ay[nbodies/2], az[nbodies/2],
	  pot[nbodies/2], potext[nbodies/2]);
*/

				/* Write output */
  out_put(n);
}
