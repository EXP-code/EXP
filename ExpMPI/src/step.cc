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

  //========================
  // Advance using leapfrog 
  //========================

				// Position step
  incr_position();

				// Compute acceleration
  comp.compute_potential();

				// Velocity step
  incr_velocity();

				// Write output
  output.Run(n);

				// Load balance
  comp.load_balance();
}
