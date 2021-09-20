                                // System libs
#include <unistd.h>

                                // C++/STL headers
#include <string>

				// Boost random generator
#include <boost/random/mersenne_twister.hpp>

				// Option parsing
#include <boost/program_options.hpp>
namespace po = boost::program_options;

                                // EXP support
#include <global.H>
#include <EmpCylSL.H>
#include <localmpi.H>

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

  std::string eof, tag;
  double rmax, zmax;
  int nout, mmax, norder;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("eof,i",      po::value<std::string>(&eof)->default_value(".eof.cache.file"),
     "the EOF cache file")
    ("out,o",      po::value<std::string>(&tag)->default_value("eof_basis"),
     "output prefix for basis functions")
    ("Rmax,R",     po::value<double>(&rmax)->default_value(0.05), 
     "Extent in cylindrical radius")
    ("Zmax,Z",     po::value<double>(&zmax)->default_value(0.005), 
     "Extent in vertical distance above and below plane")
    ("mmax,m",     po::value<int>(&mmax)->default_value(6), 
     "maximum azimuthal order")
    ("nmax,n",     po::value<int>(&norder)->default_value(18), 
     "maximum radial order")
    ("nout,N",     po::value<int>(&nout)->default_value(40), 
     "number of grid points in each dimension")
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


  std::ifstream in(eof);
  if (not in.good()) {
    if (myid==0)
      std::cerr << std::endl << argv[0] << ": error opening eof file <"
		<< eof << ">" << std::endl;

    MPI_Finalize();
    exit(-1);
  }
  in.close();

  int nmax = 64;
  int lmax = 64;
  int nord = 24;
  double acyl = 0.01;
  double hcyl = 0.001;

  EmpCylSL test(nmax, lmax, mmax, nord, acyl, hcyl);

  test.read_eof_file(eof);

  int MMAX = test.get_mmax();
  int NORD = test.get_order();

  for (int M=0; M<std::min<int>(mmax, MMAX); M++) {
    for (int N=0; N<std::min<int>(norder, NORD); N++) {
      test.dump_images_basis_pca(tag, rmax, zmax, nout, nout, M, N, 0);
    }
  }

  MPI_Finalize();
  return 0;
}

