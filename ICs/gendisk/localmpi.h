#ifndef _localmpi_h
#define _localmpi_h 1

#include <mpi.h>

void local_init_mpi(int argc, char **argv);


#ifdef IS_INIT_MPI
				/* MPI variables */
int numprocs, slaves, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

#else
				/* MPI variables */
extern int numprocs, slaves, myid, proc_namelen;
extern char* processor_name;

#endif // IS_INIT_MPI

#endif // _localmpi_h
