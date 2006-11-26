/*
  Routines to call to finish up simulation

  Global flag "finish" is set to 1
*/

#include "expand.h"
#include <global.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void clean_up(void)
{
				// Call for final output to files
  output.Run(this_step, true);

				// Debug
  cout << "Process " << myid << ": about to exit mpi . . . " << endl;

  MPI_Barrier(MPI_COMM_WORLD);
#ifdef MPE_PROFILE
  sprintf(file, "expand_mpe.%s", logfile);
  MPE_Finish_log(file);
#endif

  MPI_Finalize();

}
