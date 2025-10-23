                                // System libs
#include <unistd.h>

                                // C++/STL headers
#include <string>
                                // EXP classes
#include "libvars.H"
#include "EmpCylSL.H"
#include "localmpi.H"


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

