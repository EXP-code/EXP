/*
  Routines to call to finish up simulation

  Global flag "finish" is set to 1
*/

#include <expand.h>
#include <global.H>
#include <OutputContainer.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void clean_up(void)
{
				// Call for final output to files
  output.Run(this_step, true);

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid==0)  cerr << setfill('-') << setw(70) << "-" << endl
		     << setfill(' ')
		     << "Process " << setw(4) << right << myid 
		     << " on " << processor_name
		     << "   pid=" << getpid()
		     << "   MASTER NODE\t Exiting EXP\n";
  MPI_Barrier(MPI_COMM_WORLD);
  for (int j=1; j<numprocs; j++) {
    if (myid==j) cerr << "Process " << setw(4) << right << myid 
		      << " on " << processor_name
		      << "   pid=" << getpid()
		      << "   rank in SLAVE: " << j << "\t Exiting EXP\n";
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (myid==0)  cerr << setfill('-') << setw(70) << "-" << endl
		     << setfill(' ') << endl;

  delete parse;

  MPI_Finalize();

}
