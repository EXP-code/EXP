#define IS_INIT_MPI

#include <iostream.h>
#include <localmpi.h>

void local_init_mpi(int argc, char **argv)
{

  /*===================*/
  /* MPI preliminaries */
  /*===================*/

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  /* Debug id */
  int n;

  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_rank ( world_group, &n );
  cerr.form("Process %d on %s    rank in group: %d\n", 
	    myid, processor_name, n);

  MPI_Barrier(MPI_COMM_WORLD);
}


