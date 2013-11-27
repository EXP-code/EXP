using namespace std;

#include <unistd.h>
#include <sys/types.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <localmpi.h>

//
// MPI variables
//
MPI_Comm MPI_COMM_SLAVE;
MPI_Group world_group, slave_group;
int numprocs, slaves, myid, proc_namelen;
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
    vector<int> nslaves (slaves);

    for (int n=1; n<numprocs; n++) nslaves[n-1] = n;
    MPI_Group_incl(world_group, slaves, &nslaves[0], &slave_group);
    MPI_Comm_create(MPI_COMM_WORLD, slave_group, &MPI_COMM_SLAVE);
  }

  //=========================
  // Print node assignment
  //=========================

  if (myid==0) cerr << string(72, '-') << endl
		    << " Node, process, and communicator assignment" 
		    << endl << string(72, '-') << endl 
		    << right << setw(8)  << "Node #"
		    << " | " << setw(20) << "Hostname"
		    << " | " << setw(8)  << "PID"
		    << " | " << setw(10) << "Status"
		    << endl
		    << string(72, '-') << endl;

  for (int n=0; n<numprocs; n++) {

    if (n == myid) {

      cerr << right << setw(8)  << myid 
	   << " | " << setw(20) << processor_name 
	   << " | " << setw(8)  << getpid();
      
      if (n) {
	int m; MPI_Group_rank ( slave_group, &m );
	cerr << " | " << setw(10) << "SLAVE: " << m << endl;
      } else {
	cerr << " | " << setw(10) << "MASTER" << endl;
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
      ostringstream sout;
      sout << "mpi_debug." <<  myid;
      mpi_debug.open(sout.str().c_str());
      break;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

