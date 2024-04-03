#ifndef _localmpi_H
#define _localmpi_H

#include <mpi.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

extern MPI_Comm MPI_COMM_WORKER;
extern MPI_Group world_group, worker_group;
extern int numprocs, workers, myid, proc_namelen;
extern char processor_name[MPI_MAX_PROCESSOR_NAME];
extern std::ofstream mpi_debug;

//! Initialize MPI
void local_init_mpi(int argc, char **argv);

//! For MPI debugging
inline void mpi_report_location(const std::string& msg)
{
  for(int n=0; n<numprocs; n++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (n==myid) std::cout << "Process " << myid << ": " << msg << std::endl;
  }
}

#endif
