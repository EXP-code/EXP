                                // C++/STL headers
#include <string>

				// Option parsing
#include <cxxopts.H>

                                // MDW classes
#include <EmpCylSL.H>
#include <localmpi.H>


int 
main(int argc, char **argv)
{
  std::string eof, tag;
  double rmin, zmax;
  int nout;

 cxxopts::Options options(argv[0], "Dump the entire disk orthgonal function file in ascii");

 options.add_options()
   ("h,help", "produce this help message")
   ("i,eof", "the EOF cache file",
    cxxopts::value<std::string>(eof)->default_value(".eof.cache.file"))
   ("o,out", "output prefix for basis functions")
   cxxopts::value<std::string>(tag)->default_value("eof_basis")
   ("R,Rmax", "Extent in cylindrical radius",
    cxxopts::value<double>(rmax)->default_value("0.05"))
   ("Z,Zmax", "Extent in vertical distance above and below plane",
    cxxopts::value<double>(zmax)->default_value("0.005"))
   ("n,nout", "number of grid points in each dimension",
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
    std::cout << options.help() << "\n";
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

