/*
  Routines to call to finish up simulation

  Global flag "finish" is set to 1
*/

#include "expand.h"

void out_put(int);

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void clean_up(void)
{
				/* Set finish flag */
  finish = 1;
				/* Call for final output to files */
  out_put(nsteps+1);

  MPI_Barrier(MPI_COMM_WORLD);
#ifdef MPE_PROFILE
  sprintf(file, "expand_mpe.%s", logfile);
  MPE_Finish_log(file);
#endif
  MPI_Finalize();

}
