#include <unistd.h>		// for getpid()
#include <sys/types.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <localmpi.H>

//
// MPI variables
//
MPI_Comm MPI_COMM_SLAVE;
MPI_Group world_group, slave_group;
int numprocs=1, slaves, myid=0, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];
std::ofstream mpi_debug;

void local_init_mpi(int argc, char **argv)
{
  //===================
  // MPI preliminaries
  //===================

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  //=========================
  // Make SLAVE communicator
  //=========================

  slaves = numprocs - 1;

  if (slaves) {
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    std::vector<int> nslaves (slaves);

    for (int n=1; n<numprocs; n++) nslaves[n-1] = n;
    MPI_Group_incl(world_group, slaves, &nslaves[0], &slave_group);
    MPI_Comm_create(MPI_COMM_WORLD, slave_group, &MPI_COMM_SLAVE);
  }

  //=========================
  // Print node assignment
  //=========================

  if (myid==0)
    std::cout << std::string(80, '%') << std::endl
	      << std::setfill('%') << std::setw(80) << std::left
	      << "%%%%% Node, process, and communicator assignment " 
	      << std::endl << std::string(80, '%') << std::endl
	      << std::setfill(' ')
	      << std::right << std::setw(8)  << "Node #"
	      << " | " << std::setw(20) << "Hostname"
	      << " | " << std::setw(8)  << "PID"
	      << " | " << std::setw(10) << "Status"
	      << std::endl << std::string(80, '%') << std::endl << std::left;

  for (int n=0; n<numprocs; n++) {

    if (n == myid) {

      std::cout << std::right << std::setw(8)  << myid 
		<< " | " << std::setw(20) << processor_name 
		<< " | " << std::setw(8)  << getpid();
      
      if (n) {
	int m; MPI_Group_rank ( slave_group, &m );
	std::cout << " | " << std::setw(10) << "WORKER " << m << std::endl;
      } else {
	std::cout << " | " << std::setw(10) << "ROOT" << std::endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
    
  //=========================
  // Debugging . . . 
  //=========================

  for (int n=0; n<argc; n++) {
    std::string token(argv[n]);
    if (token.compare("--mpi-debug") == 0) {
      std::ostringstream sout;
      sout << "mpi_debug." <<  myid;
      mpi_debug.open(sout.str().c_str());
      break;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

