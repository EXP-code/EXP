#define MAIN 1
#include "expand.h"

void do_step(int);
void clean_up(void);

#include <sys/time.h>
#include <sys/resource.h>

/*
static char rcsid[] = "$Id$";
*/

void MPL_parse_args(int argc, char** argv);

int main(int argc, char** argv)
{
  const int hdbufsize=1024;
  char hdbuffer[hdbufsize];

  int *nslaves, n, retdir, retdir0;
  MPI_Group world_group, slave_group;

#ifdef DEBUG
  sleep(20);
#endif

  //===================
  // MPI preliminaries 
  //===================

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

				// Make SLAVE group 
  slaves = numprocs - 1;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  nslaves = new int [slaves];
  if (!nslaves) {
    cerr << "main: problem allocating <nslaves>\n";
    MPI_Abort(MPI_COMM_WORLD, 10);
    exit(-1);
  }
  for (n=1; n<numprocs; n++) nslaves[n-1] = n;
  MPI_Group_incl(world_group, slaves, nslaves, &slave_group);
  MPI_Comm_create(MPI_COMM_WORLD, slave_group, &MPI_COMM_SLAVE);
  delete [] nslaves;
    
  // Debug id 
  MPI_Group_rank ( slave_group, &n );
  cerr << "Process " << myid << " on " << processor_name
       << "   rank in SLAVE: " << n << "\n";
  
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef MPE_PROFILE
  MPE_Init_log();

  if (myid == 0) {
    MPE_Describe_state(1, 2, "Distribute particles", "red:dimple3" );
    MPE_Describe_state(3, 4, "Gather particles", "green:dllines3" );
    MPE_Describe_state(5, 6, "Gather coefs", "cyan:hlines2" );
    MPE_Describe_state(7, 8, "Distribute coefs", "yellow:drlines4" );
    MPE_Describe_state(9, 10, "Compute coefs", "magenta:vlines3" );
    MPE_Describe_state(11, 12, "Compute forces", "orange3:gray" );
    MPE_Describe_state(13, 14, "Advance time", "purple:boxes" );
    MPE_Describe_state(15, 16, "Send energies", "blue:dllines4" );
  }
#endif


#ifdef DEBUG
  //  sleep(20);
#endif


  //================
  // Print welcome  
  //================

  if (myid==0) {
    cout << endl << "This is " << PACKAGE << " " << VERSION
	 << " " << version_id << endl << endl;
  }


  //============================
  // Parse command line:        
  // broadcast to all processes 
  //============================

  MPL_parse_args(argc, argv);

  
  //========================
  // Change to desired home 
  // directory              
  //========================

  if (use_cwd) {
				// Get Node 0 working directory 
    if (myid == 0) getcwd(hdbuffer, (size_t)hdbufsize);
    MPI_Bcast(hdbuffer, hdbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);

    homedir.erase(homedir.begin(), homedir.end());
    homedir = hdbuffer;
    if (myid == 0) cout << "Process 0: homedir=" << homedir << "\n";
  }
  
  retdir = chdir(homedir.c_str());
  if (retdir) {
    cerr << "Process " << myid << ": could not change to home directory "
	 << homedir << "\n";
    retdir = 1;
  }
				// For exit if some nodes can't find 
				// their home 
  MPI_Allreduce(&retdir, &retdir0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (retdir0) {
    MPI_Finalize();
    exit(-1);
  }

  //=======
  // DEBUG 
  //=======
  if (myid) {
    getcwd(hdbuffer, (size_t)hdbufsize);
    cout << "Process " << myid << ": homedir=" << hdbuffer << "\n";
  }


  //================
  // Nice process ? 
  //================

  if (NICE>0) setpriority(PRIO_PROCESS, 0, NICE);


  //==============================================
  // Read in points and initialize expansion grid 
  //==============================================

  begin_run();


  //===========
  // MAIN LOOP 
  //===========

  for (this_step=1; this_step<=nsteps; this_step++) do_step(this_step);


  //===========
  // Finish up 
  //===========

  clean_up();

  return 0;
}


