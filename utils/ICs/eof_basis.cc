                                // System libs
#include <unistd.h>

                                // C++/STL headers
#include <string>

				// Option parsing
#include <cxxopts.H>

                                // EXP support
#include <global.H>
#include <EmpCylSL.H>
#include <localmpi.H>

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

  cxxopts::Options options(argv[0], "Dump the entire disk orthgonal function file in ascii");

  options.add_options()
    ("h,help", "produce this help message")
    ("i,eof", "the EOF cache file",
     cxxopts::value<std::string>(eof)->default_value(".eof.cache.file"))
    ("o,out", "output prefix for basis functions",
     cxxopts::value<std::string>(tag)->default_value("eof_basis"))
    ("R,Rmax", "Extent in cylindrical radius",
     cxxopts::value<double>(rmax)->default_value("0.05"))
    ("Z,Zmax", "Extent in vertical distance above and below plane",
     cxxopts::value<double>(zmax)->default_value("0.005"))
    ("m,mmax", "maximum azimuthal order",
     cxxopts::value<int>(mmax)->default_value("6"))
    ("n,nmax", "maximum radial order",
     cxxopts::value<int>(norder)->default_value("18"))
    ("N,nout", "number of grid points in each dimension",
     cxxopts::value<int>(nout)->default_value("40"))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0)
      std::cout << options.help() << std::endl;

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

