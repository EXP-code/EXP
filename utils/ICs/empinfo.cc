                                // System libs
#include <unistd.h>

                                // C++/STL headers
#include <string>
				// Boost random generator
#include <boost/random/mersenne_twister.hpp>

                                // MDW classes
#include <EmpCylSL.H>
#include <localmpi.H>

				// Global variables
int nthrds         = 1;
int this_step      = 0;
unsigned multistep = 0;
unsigned maxlev    = 100;
int mstep          = 1;
int Mstep          = 1;
char threading_on  = 0;
double tpos        = 0.0;
double tnow        = 0.0;

pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
string outdir, runtag;
boost::mt19937 random_gen;

void usage(char *prog)
{
  cout << setw(70) << setfill('-') << '-' << endl;
  cout << "Read and print the header file from a disk orthgonal function file" << endl;
  cout << setw(70) << setfill('-') << '-' << endl;
  cout << "Usage: " << prog << " emp_file" << endl;
  cout << setw(70) << setfill('-') << '-' << endl;
  exit(-1);
}

int 
main(int argc, char **argv)
{
  //====================
  // Inialize MPI stuff
  //====================

  local_init_mpi(argc, argv);
  
  //====================
  // Parse command line 
  //====================

  if (argc != 2) usage(argv[0]);

  EmpCylSL test;

  string eof_file(argv[1]);
  if (!test.read_eof_header(eof_file)) {
    cout << "Error reading: " << eof_file << endl;
  }

  return 0;
}

