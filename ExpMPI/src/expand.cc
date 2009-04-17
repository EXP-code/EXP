/*
static char rcsid[] = "$Id$";
*/

#include "expand.h"

void do_step(int);
void clean_up(void);

#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>

#include <fenv.h>
#include <fpetrap.h>

//===========================================
// Clean stop on a SIGTERM or SIGHUP
//===========================================

//! Abort the time stepping and checkpoint when signaled
static int stop_signal0 = 0;
void signal_handler_stop(int sig) 
{
  if (myid==0) {
    stop_signal0 = 1;
    cout << endl << "Process 0: user signaled a STOP at step=" << this_step << " . . . quitting on next step after output" << endl;
  } else {
    cout << endl << "Process " << myid << ": user signaled a STOP but only the root process can stop me . . . continuing" << endl;
  }
}

//! Dump phase space
static int dump_signal0 = 0;
void signal_handler_dump(int sig) 
{
  if (myid==0) {
    dump_signal0 = 1;
    cout << endl << "Process 0: user signaled a DUMP at step=" << this_step
	 << endl;
  } else {
    cout << endl << "Process " << myid << ": user signaled a DUMP but only the root process can do this . . . continuing" << endl;
  }
}

//! Sync argument lists on all processes
void MPL_parse_args(int argc, char** argv);

//! Print multicomputer process info
void make_node_list(int argc, char **argv)
{
  if (myid==0)  cout << setfill('-') << setw(70) << "-" << endl 
		     << setfill(' ') << endl
		     << setw(4) << left << "#" << setw(20) 
		     << "Node name" << setw(12) << "PID" 
		     << setw(40) << "Exe name" << endl
		     << setw(4) << left << "-" << setw(20) 
		     << "---------" << setw(12) << "---" 
		     << setw(40) << "--------" << endl;

  MPI_Status stat;
  unsigned nprocn = MPI_MAX_PROCESSOR_NAME, ncmd=40;
  char *procn = new char [nprocn];
  char *cmdnm = new char [ncmd];
  long pid = getpid();
  
  strncpy(procn, processor_name, nprocn);
  strncpy(cmdnm, argv[0], ncmd);

  for (int j=0; j<numprocs; j++) {
    if (myid==0) {
      if (j>0) {
	MPI_Recv(procn, nprocn, MPI_CHAR, j, 61, MPI_COMM_WORLD, &stat);
	MPI_Recv(cmdnm,   ncmd, MPI_CHAR, j, 62, MPI_COMM_WORLD, &stat);
	MPI_Recv(&pid,       1, MPI_LONG, j, 63, MPI_COMM_WORLD, &stat);
      }

      cout << setw(4) << left << j
	   << setw(20) << string(procn)
	   << setw(12) << pid 
	   << setw(ncmd) << cmdnm << endl;

    } else if (myid==j) {
      MPI_Send(procn, nprocn, MPI_CHAR, 0, 61, MPI_COMM_WORLD);
      MPI_Send(cmdnm,   ncmd, MPI_CHAR, 0, 62, MPI_COMM_WORLD);
      MPI_Send(&pid,       1, MPI_LONG, 0, 63, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  if (myid==0)  cout << setfill('-') << setw(70) << "-" << endl
		     << setfill(' ') << endl << endl;

  delete [] procn;
  delete [] cmdnm;
}

/**
   The MAIN routine
*/
int 
main(int argc, char** argv)
{
  const int hdbufsize=1024;
  char hdbuffer[hdbufsize];

#ifdef DEBUG
  set_fpu_handler();
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); 
#endif

  int *nslaves, n, retdir, retdir0;
  MPI_Group world_group, slave_group;

  //===================
  // MPI preliminaries 
  //===================

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

				// Make SLAVE group 
#ifdef SLAVE_GROUP
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
  if (myid==0)  cerr << setfill('-') << setw(70) << "-" << endl
		     << setfill(' ')
		     << "Process " << setw(4) << right << myid 
		     << " on " << processor_name
		     << "   pid=" << getpid()
		     << "   MASTER NODE\t Ready to go!\n";
#endif
#ifdef DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  for (int j=1; j<numprocs; j++) {
    if (myid==j) cerr << "Process " << setw(4) << right << myid 
		      << " on " << processor_name
		      << "   pid=" << getpid()
		      << "   rank in SLAVE: " << j << "\t Ready to go!\n";
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (myid==0)  cerr << setfill('-') << setw(70) << "-" << endl
		     << setfill(' ') << endl;
#endif

  //====================================
  // Set signal handler on HUP and TERM
  //====================================

  if (signal(SIGTERM, signal_handler_stop) == SIG_ERR) {
    cerr << endl 
	 << "Process " << myid
	 << ": Error setting signal handler [TERM]" << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
#ifdef DEBUG
  else {
    cerr << endl 
	 << "Process " << myid
	 << ": SIGTERM error handler set" << endl;
  }
#endif

  if (signal(SIGHUP, signal_handler_dump) == SIG_ERR) {
    cerr << endl 
	 << "Process " << myid
	 << ": Error setting signal handler [HUP]" << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
#ifdef DEBUG
  else {
    cerr << endl 
	 << "Process " << myid
	 << ": SIGHUP error handler set" << endl;
  }
#endif


  //================
  // Print welcome  
  //================

  if (myid==0) {
    ostringstream sout;
    cout << endl << setw(50) << setfill('%') << '%' << endl;
    sout << "%%%%% This is " << PACKAGE_STRING << " ";
    cout << left << setw(50) << sout.str() << endl;
    cout << setw(50) << "%" << setfill(' ') << endl << endl;
  }

  //====================================
  // Make node PID list
  //====================================

  make_node_list(argc, argv);

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
    if (myid == 0) cout << "main: working directory is <" << homedir << ">\n";
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
#if 0
  if (myid) {
    getcwd(hdbuffer, (size_t)hdbufsize);
    cout << "Process " << myid << ": homedir=" << hdbuffer << "\n";
  }
#endif

  //================
  // Nice process ? 
  //================

  if (NICE>0) setpriority(PRIO_PROCESS, 0, NICE);


  //==============================================
  // Read in points and initialize expansion grid 
  //==============================================

  begin_run();


  try {

    //===========
    // MAIN LOOP 
    //===========

    for (this_step=1; this_step<=nsteps; this_step++) {

      do_step(this_step);
    
      //
      // Checking for exit time
      //

      if (chktimer.done()) {
	if (myid==0) {
	  cout << "Checkpoint timer says: quit now!" << endl;
	}
      }

      //
      // Synchronize and check for signals
      //

				// Signal will only be set after the step
      dump_signal = dump_signal0;
      stop_signal = stop_signal0;
				// Reset signals
      stop_signal0 = 0;
      dump_signal0 = 0;
				// Broadcast the signal
      MPI_Bcast(&dump_signal, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(&stop_signal, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

      if (stop_signal) {
	cout << "Process " << myid << ": have stop signal\n";
	this_step++; 
	break;
      }
    
      if (dump_signal) {
	cout << "Process " << myid << ": dump signal received,"
	     << " will dump on Step " << this_step+1 << ", continuing . . .\n";
      }

    }

  }
  catch (EXPException& e) {

    cerr << "Process " << myid << ": uncaught EXP exception" << endl
	 << e.getErrorMessage() << endl;

				// Try to force all process to exit!
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  //===========
  // Finish up 
  //===========

  clean_up();

  return 0;
}


