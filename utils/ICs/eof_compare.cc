                                // System libs
#include <unistd.h>

                                // C++/STL headers
#include <string>

				// Option parsing
#include <boost/program_options.hpp>
namespace po = boost::program_options;

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
  //====================
  // Inialize MPI stuff
  //====================

  local_init_mpi(argc, argv);
  
  //====================
  // Parse command line 
  //====================

  std::string eof1, eof2;
  int nmax, lmax, mmax, norder;
  double acyl, hcyl;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help",
     "produce this help message")
    ("eof1,1",     po::value<std::string>(&eof1)->default_value(".eof.cache.file1"),
     "First EOF file")
    ("eof2,2",     po::value<std::string>(&eof2)->default_value(".eof.cache.file2"),
     "Second EOF file")
    ("mmax,m",     po::value<int>(&mmax)->default_value(64),
     "maximum spherical azimuthal order")
    ("nmax,n",     po::value<int>(&nmax)->default_value(64),
     "maximum spherical radial order")
    ("norder,N",   po::value<int>(&norder)->default_value(18), 
     "maximum cylindrical radial order")
    ("ascale,a",   po::value<double>(&acyl)->default_value(0.01), 
     "disk scale length")
    ("hscale,h",   po::value<double>(&hcyl)->default_value(0.001), 
     "disk scale height")
     ;
  
     po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0)
      std::cout << "Option error: " << e.what() << std::endl;

    MPI_Finalize();
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0)
      std::cout << desc << std::endl;

    MPI_Finalize();
    return 1;
  }


  std::ifstream in1(eof1);
  if (not in1.good()) {
    if (myid==0)
      std::cerr << std::endl << argv[0] << ": error opening eof[1] file <"
		<< eof1 << ">" << std::endl;

    MPI_Finalize();
    exit(-1);
  }
  in1.close();

  std::ifstream in2(eof2);
  if (not in2.good()) {
    if (myid==0)
      std::cerr << std::endl << argv[0] << ": error opening eof[2] file <"
		<< eof2 << ">" << std::endl;

    MPI_Finalize();
    exit(-1);
  }
  in2.close();

  EmpCylSL::RMIN        = 0.001;
  EmpCylSL::RMAX        = 20.0;
  EmpCylSL::NUMX        = 128;
  EmpCylSL::NUMY        = 64;
  EmpCylSL::CMAP        = true;
  EmpCylSL::DENS        = true;
  EmpCylSL::VFLAG       = 26;

  EmpCylSL::CACHEFILE = eof1;

  EmpCylSL test1(nmax, lmax, mmax, norder, acyl, hcyl);

  bool cache_ok = test1.read_cache();

  if (!cache_ok) {
    if (myid==0) {		// Diagnostic output . . .
      std::cerr << "Can not read explicitly specified EOF file <"
		<< EmpCylSL::CACHEFILE << ">" << std::endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  EmpCylSL::CACHEFILE = eof2;

  EmpCylSL test2(nmax, lmax, mmax, norder, acyl, hcyl);

  cache_ok = test1.read_cache();

  if (!cache_ok) {
    if (myid==0) {		// Diagnostic output . . .
      std::cerr << "Can not read explicitly specified EOF file <"
		<< EmpCylSL::CACHEFILE << ">" << std::endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  if (myid==0) {
    test1.compare_basis(&test2);
  }

  MPI_Finalize();
  return 0;
}

