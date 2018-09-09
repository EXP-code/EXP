#include "expand.h"

void do_step (int);
void clean_up(void);

#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>

#include <fenv.h>
#include <fpetrap.h>

#ifdef USE_GPTL
#include <gptl.h>
#endif

#include <BarrierWrapper.H>
#include <FileUtils.H>


//===========================================
// Handlers defined in exputil/stack.cc
//===========================================

extern void mpi_print_trace(const string& routine, const string& msg,
			    const char *file, int line);

extern void mpi_gdb_print_trace(int sig);

extern void mpi_gdb_wait_trace(int sig);

//===========================================
// A signal handler to produce a traceback
//===========================================

void set_fpu_trace_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);
}

//===========================================
// A signal handler to produce stop and wait
//===========================================

void set_fpu_gdb_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_wait_trace);
}


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

  // Check for barrier wrapper failure to provide some additional
  // debugging info to help with interpreting the backtrace . . .
  //
  ostringstream sout;
  if (BarrierWrapper::inOper) {
    sout << "BarrierWrapper label is:" 
	 << left << setw(60) << BarrierWrapper::lbOper;
    if (BarrierWrapper::flOper.size())
      sout << " ==> called from " << BarrierWrapper::flOper << ":" 
	   << BarrierWrapper::lnOper;
  } else
    sout << "process called abort";

  // Print the backtrace to per process files or stderr
  //
  mpi_print_trace("signal_handler_stop", sout.str(), __FILE__, __LINE__);
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
		     << setw(40) << "Executable" << endl
		     << setw(4) << left << "-" << setw(20) 
		     << "---------" << setw(12) << "---" 
		     << setw(40) << "----------" << endl;

  MPI_Status stat;
  unsigned nprocn = MPI_MAX_PROCESSOR_NAME, ncmd=40;
  char *procn = new char [nprocn];
  char *cmdnm = new char [ncmd];
  long pid = getpid();
  
  strncpy(procn, processor_name, nprocn);
  strncpy(cmdnm, argv[0], ncmd);

  if (myid==0) {
    for (int j=0; j<numprocs; j++) {
      if (j) {
	MPI_Recv(procn, nprocn, MPI_CHAR, j, 61, MPI_COMM_WORLD, &stat);
	MPI_Recv(cmdnm,   ncmd, MPI_CHAR, j, 62, MPI_COMM_WORLD, &stat);
	MPI_Recv(&pid,       1, MPI_LONG, j, 63, MPI_COMM_WORLD, &stat);
      }
      cout << setw(4)  << left << j
	   << setw(20) << string(procn)
	   << setw(12) << pid 
	   << setw(ncmd) << cmdnm << endl;
    }
  } else {
    MPI_Send(procn, nprocn, MPI_CHAR, 0, 61, MPI_COMM_WORLD);
    MPI_Send(cmdnm,   ncmd, MPI_CHAR, 0, 62, MPI_COMM_WORLD);
    MPI_Send(&pid,       1, MPI_LONG, 0, 63, MPI_COMM_WORLD);
  }
  
  if (myid==0)  cout << setfill('-') << setw(70) << "-" << endl
		     << setfill(' ') << endl << endl;

  // Make MPI datatype

#ifdef INT128
  int ityp = 2;
#else
  int ityp = 1;
#endif
  MPI_Aint kdsp = 0;
  MPI_Datatype ktyp = MPI_UNSIGNED_LONG;
  MPI_Type_create_struct(1, &ityp, &kdsp, &ktyp, &MPI_EXP_KEYTYPE);
  MPI_Type_commit(&MPI_EXP_KEYTYPE);

  // Generate node list for this node

  int  *total_ranks = new int  [sizeof(int) * numprocs];
  char *total_names = new char [MPI_MAX_PROCESSOR_NAME * numprocs];

				// Send all ranks and host names
  MPI_Allgather(&myid, 1, MPI_INT, total_ranks, 1, MPI_INT,
		MPI_COMM_WORLD);

  MPI_Allgather(processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
		total_names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
		MPI_COMM_WORLD);

				// Catagorize all ranks by host names
  for (int i=0; i<numprocs; i++) {
    char tmp[MPI_MAX_PROCESSOR_NAME];
    strncpy(tmp, &total_names[MPI_MAX_PROCESSOR_NAME*i], MPI_MAX_PROCESSOR_NAME);
    nameMap[std::string(tmp)].push_back(total_ranks[i]);
  }

				// Sort rank list (not really necessariy)
  for (auto & v : nameMap) {
    std::sort(v.second.begin(), v.second.end());
  }

				// Sibling ranks on this processor
  siblingList = nameMap[std::string(processor_name)];

  delete [] total_names;
  delete [] total_ranks;

  delete [] procn;
  delete [] cmdnm;
}

std::map<std::string, std::vector<int> >
generateNodeList()
{
  

  return nameMap;
}


int set_memlock_limits()
{
  if (rlimit_val==0) return 0;

  struct rlimit rlim;
  
  if (rlimit_val<0) {
    rlim.rlim_cur = RLIM_INFINITY;
    rlim.rlim_max = RLIM_INFINITY;
  } else {
    const rlim_t GB = 1024*1024*1024;
    rlim.rlim_cur = static_cast<rlim_t>(rlimit_val) * GB;
    rlim.rlim_max = static_cast<rlim_t>(rlimit_val) * GB;
  }
  
  return setrlimit(RLIMIT_MEMLOCK, &rlim);
}


void report_memlock_limits()
{
  struct rlimit rlim;
  ostringstream sout;

  sout << "Node " << processor_name << ", pid=" << getpid();
  if (getrlimit(RLIMIT_MEMLOCK, &rlim) == 0) {
    sout << ", using MEMLOCK limits (soft, hard) = (" 
	 << rlim.rlim_cur << ", " << rlim.rlim_max << ")";
  } else {
    sout << ", could not get RLIMIT_MEMLOCK!";
  }

  cout << sout.str() << endl << flush;
}


/**
   The MAIN routine
*/
int 
main(int argc, char** argv)
{
  const int hdbufsize=1024;
  char hdbuffer[hdbufsize];
  bool final_cmd = false;

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


  //============
  // Start GPTL
  //============

#ifdef USE_GPTL
  {
    int ret;

    ret = GPTLsetoption (GPTLoverhead,       1);
    ret = GPTLsetoption (GPTLpercent,        1);
    ret = GPTLsetoption (GPTLabort_on_error, 0);

    ret = GPTLinitialize();
    ret = GPTLstart("main");
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

  //============================
  // Trap floating point errors
  // by installing user handler
  //============================

  if (fpe_trap ) set_fpu_handler();
  if (fpe_trace) set_fpu_trace_handler();
  if (fpe_wait ) set_fpu_gdb_handler();

  //========================
  // Change to desired home 
  // directory              
  //========================

  if (use_cwd) {
				// Get Node 0 working directory 
    if (myid == 0) auto ret = getcwd(hdbuffer, (size_t)hdbufsize);
    MPI_Bcast(hdbuffer, hdbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);

    homedir.clear();
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

  //==================
  // Barrier debugging
  //==================

  barrier = new BarrierWrapper(MPI_COMM_WORLD, barrier_label);
  if (barrier_check) barrier->on();
  else               barrier->off();
  if (barrier_light) barrier->setLightWeight();
  else               barrier->setHeavyWeight();
  if (barrier_quiet) BarrierWrapper::verbose       = false;
  else               BarrierWrapper::verbose       = true;
  if (barrier_extra) BarrierWrapper::extra_verbose = true;
  else               BarrierWrapper::extra_verbose = false;
  if (barrier_debug) BarrierWrapper::debugging     = true;
  else               BarrierWrapper::debugging     = false;

  //================
  // Nice process ? 
  //================

  if (NICE>0) setpriority(PRIO_PROCESS, 0, NICE);


  //===================
  // Set memory limits
  //===================

  if (set_memlock_limits()) report_memlock_limits();

  //==============================================
  // Sleep loop for debugging
  //==============================================

  if (!main_wait or myid==0) {
    while (debug_wait) {
      sleep(5);
    }
  }

  //==============================================
  // Read in points and initialize expansion grid 
  //==============================================

  (*barrier)("Expand: BEFORE begin_run", __FILE__, __LINE__);
  begin_run();
  (*barrier)("Expand: AFTER begin_run", __FILE__, __LINE__);

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
	  cout << "Restart command is: " << restart_cmd << endl;
	  quit_signal = 1;
	}
      }
      MPI_Bcast(&quit_signal, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);      
      if (quit_signal)  break;

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

    std::cerr << "Process " << myid << ": uncaught EXP exception" << std::endl
	      << e.getErrorMessage() << std::endl;

				// Try to force all process to exit!
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  //=============
  // Finish GPTL
  //=============

#ifdef USE_GPTL
  {
    int ret, cnt = 0;
    const int safety = 10000;
    ret = GPTLstop("main");
    ostringstream sout;
    sout << runtag << "_timing." << myid;
    string tfile = sout.str();
    if (FileExists(tfile)) {
      do {
	sout.str("");
	sout << tfile << ".bak." << cnt++;
      } while (FileExists(sout.str().c_str()) && cnt<safety);
      FileRename(tfile, sout.str().c_str());
    }
    ret = GPTLpr_file(tfile.c_str());
  }
#endif

  //=================
  // A clean shutdown
  //=================

  clean_up();

  //=================
  // Epilogue command
  //=================

  if (quit_signal && myid==0) {
    std::cout << "Executing the epilogue command: " << restart_cmd << endl;
    if (system(restart_cmd.c_str()) == -1) {
      std::cerr << "In MAIN, error executing the restart command <" 
		<< restart_cmd << ">" << std::endl;
    }
  }

  return 0;
}


