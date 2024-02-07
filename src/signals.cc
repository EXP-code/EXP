#include <signal.h>
#include <unistd.h>

#include <fpetrap.h>
#include "expand.H"

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
