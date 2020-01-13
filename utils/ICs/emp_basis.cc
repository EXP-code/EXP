                                // System libs
#include <unistd.h>

                                // C++/STL headers
#include <string>

				// Option parsing
#include <boost/program_options.hpp>

                                // MDW classes
#include <EmpCylSL.h>
#include <localmpi.h>

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
  std::string eof, tag;
  double rmin, zmax;
  int nout;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("eof,i",      po::value<std::string>(&eof)->default_value(".eof.cache.file")
     "the EOF cache file")
    ("out,o",      po::value<std::string>(&tag)->default_value("eof_basis")
     "output prefix for basis functions")
    ("Rmax,R",     po::value<double>(&rmax)->default_value(0.05), 
     "Extent in cylindrical radius")
    ("Zmax,Z",     po::value<double>(&zmax)->default_value(0.005), 
     "Extent in vertical distance above and below plane")
    ("nout,n",     po::value<int>(&nout)->default_value(40), 
     "number of grid points in each dimension")
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }


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

