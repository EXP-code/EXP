using namespace std;

#include <iostream>
#include <localmpi.h>

// MPI variables
//
MPI_Comm MPI_COMM_SLAVE;
int numprocs, slaves, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

void local_init_mpi(int argc, char **argv)
{
  MPI_Group world_group, slave_group;
  int n;

  //===================
  // MPI preliminaries
  //===================

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  // Make SLAVE communicator
  //
  slaves = numprocs - 1;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  int *nslaves = new int [slaves];

  for (n=1; n<numprocs; n++) nslaves[n-1] = n;
  MPI_Group_incl(world_group, slaves, nslaves, &slave_group);
  MPI_Comm_create(MPI_COMM_WORLD, slave_group, &MPI_COMM_SLAVE);
  delete [] nslaves;
    
  // Debug id
  //
  MPI_Group_rank ( slave_group, &n );
  cerr << "Process " << myid << " on " << processor_name 
       << "   rank in SLAVE: "
       << n << endl;

  MPI_Barrier(MPI_COMM_WORLD);
}
