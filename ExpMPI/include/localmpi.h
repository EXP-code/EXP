#include <mpi.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif


#ifdef IS_MAIN
				/* MPI variables */
int numprocs, slaves, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

#else
				/* MPI variables */

extern int numprocs, slaves, myid, proc_namelen;
extern char processor_name[MPI_MAX_PROCESSOR_NAME];

#endif /* IS_MAIN */

