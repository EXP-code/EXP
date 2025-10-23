#include <unistd.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <signal.h>
#include <fenv.h>

// For debugging . . . unwinds stack and writes to output stream with
// symbols if available

#include <mpi.h>

#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <cstdlib>
#include <stdexcept>

#include "libvars.H"
using namespace __EXP__;

//! Traceback type
bool gdb_trace     = true;


//! Use native Linux backtrace functions to get stack info
void print_trace(std::ostream& out, const char *file, int line)
{
  const size_t max_depth = 100;
  size_t       stack_depth;
  void        *stack_addrs[max_depth];
  char       **stack_strings;

  //
  // These are GNU C extensions and therefore are not going to be
  // portable, alas.
  //
  stack_depth   = backtrace        (stack_addrs, max_depth  );
  stack_strings = backtrace_symbols(stack_addrs, stack_depth);
  
  out << std::setfill('-') << std::setw(80) << '-' << std::endl;
  out << std::setfill(' ');

  if (file) out << "Call stack from " << file << ":" << line << std::endl;
  
  if (0) {
    for (size_t i = 1; i < stack_depth; i++) {
      out << "    " << stack_strings[i] << std::endl;
    }
  }
  
  for (size_t i = 1; i < stack_depth; i++) {
    //
    // 4 x 80 character lines worth.  I suppose some template names
    // may be larger . . .
    //
    size_t sz = 320;
				// We need to use malloc/free for
				// these functions, sigh . . .
    char *function = static_cast<char *>(malloc(sz));
    char *begin = 0, *offset = 0, *end = 0;
    //
    // Find the parentheses and address offset surrounding the mangled
    // name
    //
    for (char *j = stack_strings[i]; *j; ++j) {
      if (*j == '(')
	begin = j;
      else if (*j == '+')
	offset = j;
      else if (*j == ')' && offset) {
	end = j;
      }
    }
    if (begin && offset && end && begin < offset) {
      *begin++  = '\0';
      *offset++ = '\0';
      *end      = '\0';

      std::ostringstream sout;

      //
      // Found the mangled name, now in [begin, end)
      //	
      int status;
      char *ret = abi::__cxa_demangle(begin, function, &sz, &status);
      if (ret) {
	//
	// Return value may be a realloc() of the input
	//
	function = ret;

	sout << "  " << stack_strings[i] << " : "
	     << ret << "()+" << offset << std::endl;
      }
      else {
	//
	// Demangling failed, format it as a C function with no args
	//

	sout << "  " << stack_strings[i] << " : "
	     << begin << "()+" << offset << std::endl;
      }
      out << sout.str();
    } else {
				// Didn't find the mangled name, just
				// print the whole line
      out << "    " << stack_strings[i] << std::endl;
    }
    free(function);		// malloc()ed above
  }

  free(stack_strings);		// malloc()ed by backtrace_symbols
  
  out << std::setfill('-') << std::setw(80) << '-' << std::endl
      << std::setfill(' ') << std::flush;
}


void mpi_print_trace(const std::string& routine, const std::string& msg,
		     const char *file, int line)
{
  //
  // Look for active MPI environment
  //
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  std::cerr << routine;
  if (numprocs>1) std::cerr << " [mpi_id=" << myid << "]";
  std::cerr << ": "<< msg << std::endl;
  
  std::ostringstream ostr;
  ostr << outdir << runtag << "." << "traceback.";
  if (numprocs>1) ostr << myid;
  else            ostr << "info";
  
  std::ofstream tb(ostr.str().c_str());

  //
  // Print out all the frames
  //
  if (tb.good()) {
    std::cerr << routine;
    if (numprocs>1) std::cerr << " [mpi_id=" << myid << "]";
    std::cerr << ": see <" << ostr.str() << "> for more info" << std::endl;
    
    print_trace(tb,        0, 0);
  } else {
    print_trace(std::cerr, 0, 0);
  }
}


std::string exec(const std::string& cmd)
{
  std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);

  if (!pipe) throw std::runtime_error("popen() failed!");

  char buffer[128];
  std::string result = "";
  while (!feof(pipe.get())) {
    if (fgets(buffer, 128, pipe.get()) != NULL)
      result += buffer;
  }
  return result;
}

void print_gdb_backtrace(std::ostream & out)
{
  pid_t pid = getpid();
  char name[512];
  ssize_t nlen = readlink("/proc/self/exe", name, 511);
  
  std::ostringstream command;
  command << "gdb --batch -nx -ex thread -ex bt -p " << pid << " " <<  name;
  
  std::string result = exec(command.str());
  
  out << result;
}

void mpi_gdb_print_trace(int sig)
{
  //
  // Look for active MPI environment
  //
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  std::ostringstream sout;
  sout << "FPE error, node=" << myid << ": signal " << sig;

  std::cerr << sout.str();
  if (numprocs>1) std::cerr << " [mpi_id=" << myid << "]";
  std::cerr << std::endl;
  
  std::ostringstream ostr;
  ostr << outdir << runtag << "." << "traceback.";
  if (numprocs>1) ostr << myid;
  else            ostr << "info";
  
  std::ofstream tb(ostr.str().c_str());

  //
  // Print out all the frames
  //
  if (tb.good()) {
    if (numprocs>1) std::cerr << " [mpi_id=" << myid << "]";
    std::cerr << ": see <" << ostr.str() << "> for more info" << std::endl;
    tb << sout.str() << std::endl;
    if (gdb_trace)
      print_gdb_backtrace(tb       );
    else
      print_trace(tb, 0, 0);
  } else {
    if (gdb_trace)
      print_gdb_backtrace(std::cerr);
    else
      print_trace(std::cerr, 0, 0);
  }

  tb.close();

  exit(1);
}

void mpi_gdb_wait_trace(int sig)
{
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
  sout << "FPE error, name=" << hostname
       << ", pid=" << pid
       << ", node=" << myid << ": signal " << sig;

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
      std::cerr << hostname << "[pid=" << pid << "] waiting "
		<< sofar/60 << " minutes" << std::endl;
				// Notify every ten minutes
    } else if (sofar < 3600 and sofar % 600==0) {
      std::cerr << hostname << "[pid=" << pid << "] waiting "
		<< sofar/60 << " minutes" << std::endl;
				// Notify every thirty minutes
    } else if (sofar % 1800==0) {
      std::cerr << hostname << "[pid=" << pid << "] waiting "
		<< sofar/60 << " minutes" << std::endl;
    }      
  }
}
