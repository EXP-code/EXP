/*
  Get the whole thing going by:
  -----------------------------
  1) reading in initial phase space
  2) computing potential initially
  3) dumping first phase-space and log file
  4) write parameter file
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void out_put(int);
void adjust_scale(int);
void write_parm(void);

void begin_run(void)
{

  /*================================*/
  /* Initialize phase-space vectors */
  /*================================*/

  do_bodies_init();


  /*====================================*/
  /* Distribute particles to processors */
  /*====================================*/

  read_bodies_and_distribute();

  /*==================================*/
  /* put COM at origin in phase space */
  /*==================================*/

  if (fixpos==1) fix_positions();
  if (fixpos==2) fix_positions_by_component();


  /*===============================*/
  /* Compute initial accereration  */
  /*===============================*/

  compute_potential();


  /*================================*/
  /* Add fictitious accereration to */
  /* keep COM at origin             */
  /*================================*/

  if (fixacc) fix_acceleration();

  /*======================================*/
  /* Initial dump file (zero forces dump) */
  /*======================================*/

  out_put(0);


  /*==========================*/
  /* Find new scale parameter */
  /*==========================*/

  if (adjust) adjust_scale(0);

  /*==============================*/
  /* Increment velocity of bodies */
  /* to initialize leap-frog      */
  /*==============================*/

  init_velocity();

  /*======================*/
  /* Write parameter file */
  /*======================*/

  write_parm();
}
