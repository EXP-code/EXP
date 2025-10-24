                                // System libs
#include <unistd.h>

                                // C++/STL headers
#include <string>
#include <random>

                                // EXP library support
#include "libvars.H"
#include "EmpCylSL.H"
#include "localmpi.H"

void usage(char *prog)
{
  cout << setw(70) << setfill('-') << '-' << endl;
  cout << "Dump the entire disk orthgonal function file in ascii" << endl;
  cout << setw(70) << setfill('-') << '-' << endl;
  cout << "Usage: " << prog << " emp_file dump_file" << endl;
  cout << setw(70) << setfill('-') << '-' << endl;

  MPI_Finalize();
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

  if (argc != 3)         usage(argv[0]);
  if (argv[1][0] == '-') usage(argv[0]);

  std::ifstream in(argv[1]);
  if (not in.good()) {
    std::cerr << std::endl << argv[0] << ": error opening eof file <"
	      << argv[1] << ">" << std::endl;
    MPI_Finalize();
    exit(-1);
  }
  in.close();

  EmpCylSL test;

  string eof_file(argv[1]);
  string dmp_file(argv[2]);

  test.dump_eof_file(eof_file, dmp_file);

  MPI_Finalize();

  return 0;
}

