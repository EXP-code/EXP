#include <iostream>
#include "localmpi.h"

void local_init_mpi(int argc, char **argv)
{

  /*===================*/
  /* MPI preliminaries */
  /*===================*/

  processor_name = new char [MPI_MAX_PROCESSOR_NAME];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  /* Debug id */
  int n;

  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_rank ( world_group, &n );
  std::cerr << "Process" << myid << "on" << processor_name << "rank in group:" << n << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
}


