#ifndef _localmpi_h
#define _localmpi_h

#include <mpi.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

extern int numprocs, slaves, myid, proc_namelen;
extern char processor_name[MPI_MAX_PROCESSOR_NAME];

#include <iostream>
#include <iomanip>
#include <string>

//! For MPI debugging
inline void mpi_report_location(const std::string& msg)
{
  for(int n=0; n<numprocs; n++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (n==myid) std::cout << "Process " << myid << ": " << msg << std::endl;
  }
}

#endif
