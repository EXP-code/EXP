#ifndef _localmpi_h
#define _localmpi_h

#include <mpi.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

// We don't want to redefine these if they have already been defined
// elsewhere . . . 
//
#ifndef _global_H
extern MPI_Comm MPI_COMM_SLAVE;
extern MPI_Group world_group, slave_group;
extern int numprocs, slaves, myid, proc_namelen;
extern char processor_name[MPI_MAX_PROCESSOR_NAME];
extern std::ofstream mpi_debug;
#endif

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
