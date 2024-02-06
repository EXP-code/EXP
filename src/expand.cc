#include "expand.H"

void do_step (int);
void clean_up(void);

#include <signal.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>

#include <fpetrap.h>

#ifdef USE_GPTL
#include <gptl.h>
#endif

#ifdef HAVE_LIBSLURM
#include <slurm/slurm.h>
#endif

#include <BarrierWrapper.H>
#include <FileUtils.H>


//===========================================
// Helper defined in parse.cc
//===========================================

extern void exp_version();

//===========================================
// Handlers defined in exputil/stack.cc
//===========================================

extern void print_trace(std::ostream& out,
			const char *file, int line);

extern void mpi_print_trace(const string& routine, const string& msg,
			    const char *file, int line);

extern void mpi_gdb_print_trace(int sig);

extern void mpi_gdb_wait_trace(int sig);

//===========================================
// A signal handler to trap invalid FP only
//===========================================

void set_fpu_invalid_handler(void)
{
  // These calls are provided by glibc.  The key function
  // 'feenableexcept(/*flags*/);' has the system raise a signal that
  // may be trapped by gdb for debugging.
  // 
  // Please contribute solutions
  // for other systems and unsupported architectures if possible...

#ifdef HAVE_FE_ENABLE
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
    std::cout << "---- Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);

#endif
}

//===========================================
// A signal handler to produce a traceback
//===========================================

void set_fpu_trace_handler(void)
{
#ifdef HAVE_FE_ENABLE
  
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

#endif
}

//===========================================
// A signal handler to produce stop and wait
//===========================================

void set_fpu_gdb_handler(void)
{
#ifdef HAVE_FE_ENABLE

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
#endif
}



void exp_mpi_error_handler(MPI_Comm *communicator, int *error_code, ...)
{
  char error_string[MPI_MAX_ERROR_STRING];
  int error_string_length;

  std::cout << "exp_mpi_error_handler: entry" << std::endl;
  std::cout << "exp_mpi_error_handler: error_code = " << *error_code
	    << std::endl;
  MPI_Error_string(*error_code, error_string, &error_string_length);
  error_string[error_string_length] = '\0';

  //
  // Look for active MPI environment
  //
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  //
  // Get host name and pid
  //
  const size_t maxlen = 128;
  char hostname[maxlen];
  int hsiz = gethostname(hostname, maxlen);
  pid_t pid = getpid();

  std::ostringstream sout;
  sout << "MPI error, name=" << hostname
       << ", pid=" << pid
       << ", node=" << myid << ": signal " << error_code
       << ", error=" << error_string;

  std::cerr << sout.str();
  if (numprocs>1) std::cerr << " [mpi_id=" << myid << "]";
  std::cerr << std::endl;

  //
  // Enter sleep loop
  //
  const unsigned int waittime = 2;
  unsigned int sofar = 0;
  bool go = true;
  while (go) {
    sleep(waittime);
    sofar += waittime;
				// Notify every minute for 10 minutes
    if (sofar < 600 and sofar % 60==0) {
      std::cerr << hostname << "[pid=" << pid << "] waiting [MPI] "
		<< sofar/60 << " minutes" << std::endl;
				// Notify every ten minutes
    } else if (sofar < 3600 and sofar % 600==0) {
      std::cerr << hostname << "[pid=" << pid << "] waiting [MPI] "
		<< sofar/60 << " minutes" << std::endl;
				// Notify every thirty minutes
    } else if (sofar % 1800==0) {
      std::cerr << hostname << "[pid=" << pid << "] waiting [MPI] "
		<< sofar/60 << " minutes" << std::endl;
    }      
  }

}


//===========================================
// Clean stop on a SIGTERM or SIGHUP
//===========================================

//! Abort the time stepping and checkpoint when signaled
static int stop_signal0 = 0;

//! Dump phase space
static int dump_signal0 = 0;

void signal_handler_stop(int sig) 
{
  if (myid==0) {
    stop_signal0 = 1;
    dump_signal0 = 1;
    std::cout << std::endl
	      << "Process 0: user signaled a STOP at step=" << this_step
	      << " . . . quitting on next step after output" << std::endl;
  } else {
    std::cout << std::endl
	      << "Process " << myid
	      << ": user signaled a STOP but only the root process can stop me . . . continuing" << std::endl;
  }

  // Check for barrier wrapper failure to provide some additional
  // debugging info to help with interpreting the backtrace . . .
  //
  std::ostringstream sout;
  if (BarrierWrapper::inOper) {
    sout << "BarrierWrapper label is:" 
	 << std::left << setw(60) << BarrierWrapper::lbOper;
    if (BarrierWrapper::flOper.size())
      sout << " ==> called from " << BarrierWrapper::flOper << ":" 
	   << BarrierWrapper::lnOper;
  } else
    sout << "process called abort";

  // Print the backtrace to per process files or stderr
  //
  mpi_print_trace("signal_handler_stop", sout.str(), __FILE__, __LINE__);
}

void signal_handler_dump(int sig) 
{
  if (myid==0) {
    dump_signal0 = 1;
    std::cout << std::endl
	      << "Process 0: user signaled a DUMP at step=" << this_step
	      << std::endl;
  } else {
    std::cout << std::endl
	      << "Process " << myid << ": user signaled a DUMP but only the root process can do this . . . continuing"
	      << std::endl;
  }
}

//! Sync argument lists on all processes
void YAML_parse_args(int argc, char** argv);

//! Print multicomputer process info
void make_node_list(int argc, char **argv)
{
  if (myid==0)
    std:: cout << std::setfill('%') << std::setw(80) << "%" << std::endl 
	       << std::setfill(' ') << std::endl
	       << std::setw(4) << std::left << "#" <<std::setw(20) 
	       << "Node name" << std::setw(12) << "PID" 
	       << std::setw(40) << "Executable" << endl
	       << std::setw(4) << std::left << "-" << std::setw(20) 
	       << "---------" << std::setw(12) << "---" 
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
      std::cout << std::setw(4)  << std::left << j
		<< std::setw(20) << std::string(procn)
		<< std::setw(12) << pid 
		<< std::setw(ncmd) << cmdnm << std::endl;
    }
  } else {
    MPI_Send(procn, nprocn, MPI_CHAR, 0, 61, MPI_COMM_WORLD);
    MPI_Send(cmdnm,   ncmd, MPI_CHAR, 0, 62, MPI_COMM_WORLD);
    MPI_Send(&pid,       1, MPI_LONG, 0, 63, MPI_COMM_WORLD);
  }
  
  if (myid==0) std::cout << std::setfill('%') << std::setw(80) << "%" << std::endl
			 << std::setfill(' ') << std::endl << std::endl;

  // Make MPI datatype

#ifdef I64
  MPI_EXP_KEYTYPE = MPI_UNSIGNED_LONG_LONG;
#else
  MPI_EXP_KEYTYPE = MPI_UNSIGNED_LONG;
#endif

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

  std::cout << sout.str() << std::endl << std::flush;
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
  MPI_Errhandler errhandler;

  //===================
  // MPI preliminaries 
  //===================

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

				// Installation of MPI error handler
  if (mpi_wait) {
    MPI_Comm_create_errhandler(&exp_mpi_error_handler, &errhandler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
  }

				// Make SLAVE group 
#ifdef SLAVE_GROUP
  slaves = numprocs - 1;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  nslaves = new int [slaves];
  if (!nslaves) {
    cerr << "main: problem allocating <nslaves>\n";
    MPI_Finalize();
    exit(10);
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
    MPI_Finalize();
    exit(11);
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
    MPI_Finalize();
    exit(12);
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

  //====================================
  // Make node PID list
  //====================================

  make_node_list(argc, argv);

  //================
  // Print welcome  
  //================

  if (myid==0) exp_version();

  
  //==============================================
  // Begin exception handling
  //==============================================

  try {

    //============================
    // Parse command line:        
    // broadcast to all processes 
    //============================

    YAML_parse_args(argc, argv);

    //============================
    // Trap floating point errors
    // by installing user handler
    //============================

    if (fpe_trap ) set_fpu_invalid_handler();
    if (fpe_trace) set_fpu_trace_handler();
    if (fpe_wait ) set_fpu_gdb_handler();

    //========================
    // Change to desired home 
    // directory              
    //========================
    
    if (use_cwd) {
      // Get Node 0 working directory and broadcast to everyone
      //
      if (myid == 0) auto ret = getcwd(hdbuffer, (size_t)hdbufsize);
      MPI_Bcast(hdbuffer, hdbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);

      homedir.clear();
      homedir = hdbuffer;
      if (myid == 0)
	std::cout << "main: working directory is <" << homedir << ">"
		  << std::endl;
    }
  
    retdir = chdir(homedir.c_str());
    if (retdir) {
      std::cerr << "Process " << myid << ": could not change to home directory "
		<< homedir << std::endl;
      retdir = 1;
    }

    // For exit if some nodes can't find their home
    //
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
      std::cout << "Process " << myid << ": homedir=" << hdbuffer << std::endl;
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
    
    //=====================
    // Initialize slurm api
    //=====================
#ifdef HAVE_LIBSLURM
#if SLURM_VERSION_NUMBER > SLURM_VERSION_NUM(21,10,0)
    slurm_init(0);
#endif
#endif

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
	  std::cout << "Checkpoint timer says: quit now!" << std::endl;
	  std::cout << "Restart command is: " << restart_cmd << std::endl;
	  quit_signal = 1;
	}
      }
      MPI_Bcast(&quit_signal, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);      

      //
      // Synchronize and check for signals
      //

      // Signal will only be set after the step
      //
      dump_signal = dump_signal0;
      stop_signal = stop_signal0;

      // Reset static signals
      //
      stop_signal0 = 0;
      dump_signal0 = 0;

      // Broadcast the signals
      //
      MPI_Bcast(&dump_signal, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
      MPI_Bcast(&stop_signal, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

      // Break the step loop if quit flag is set
      //
      if (quit_signal) break;

      // Be chatty about OS signals
      //
      if (stop_signal) {
	std::cout << "Process " << myid << ": have stop signal" << std::endl;
	this_step++; 
	break;
      }
    
      if (dump_signal) {
	std::cout << "Process " << myid << ": dump signal received,"
		  << " will dump on Step " << this_step+1
		  << ", continuing . . ." << std::endl;
      }
    }
    
  }
  catch (EXPException& e) {

    if (e.getDeadlock())
      std::cerr << std::string(72, '-') << std::endl
		<< "Process " << myid << ": EXP asynchronous exception"
		<< std::endl
		<< std::string(72, '-') << std::endl
		<< e.getErrorMessage()  << std::endl
		<< std::string(72, '-') << std::endl;
    else if (myid==0)
      std::cerr << std::string(72, '-') << std::endl
		<< "EXP synchronous exception" << std::endl
		<< std::string(72, '-') << std::endl
		<< e.getErrorMessage()  << std::endl
		<< std::string(72, '-') << std::endl;

    if (traceback and (e.getDeadlock() or myid==0)) {
      print_trace(std::cerr, 0, 0);
      sleep(5);
      std::cerr << std::flush;
    }

    // Try to force all process to exit!
    if (e.getDeadlock()) MPI_Abort(MPI_COMM_WORLD, e.getErrorcode());
    else MPI_Finalize();

    exit(0);
  }
  catch (std::runtime_error& e) {
    std::cerr << "Process " << myid << ": std exception" << std::endl
	      << e.what() << std::endl;
    if (VERBOSE>4) print_trace(std::cerr, 0, 0);
    sleep(5);
    std::cerr << std::flush;

    // Try to force all process to exit!
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  catch (std::string& msg) {
    std::cerr << "Process " << myid << ": str exception" << std::endl
	      << msg << std::endl;
    sleep(5);
    if (VERBOSE>4) print_trace(std::cerr, 0, 0);
    std::cerr << std::flush;

    // Try to force all process to exit!
    MPI_Abort(MPI_COMM_WORLD, 0);
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

  //======================
  // Clean up slurm config
  //======================

#ifdef HAVE_LIBSLURM
#if SLURM_VERSION_NUMBER > SLURM_VERSION_NUM(21,10,0)
  slurm_fini();
#endif
#endif

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


