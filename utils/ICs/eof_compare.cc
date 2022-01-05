                                // C++/STL headers
#include <string>

                                // EXP support
#include <global.H>
#include <EmpCylSL.H>
#include <localmpi.H>
#include <cxxopts.H>


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
  int nmax, lmax, mmax, norder, nodd;
  double acyl, hcyl;

  cxxopts::Options options(argv[0], "Compute two orthogonal function files");
  
  options.add_options()
   ("help", "produce this help message")
   ("1,eof1", "First EOF file",
     cxxopts::value<std::string>(eof1)->default_value(".eof.cache.file1"))
   ("2,eof2", "Second EOF file",
     cxxopts::value<std::string>(eof2)->default_value(".eof.cache.file2"))
   ("l,lmax", "maximum spherical azimuthal order",
     cxxopts::value<int>(lmax)->default_value("32"))
   ("n,nmax", "maximum spherical radial order",
     cxxopts::value<int>(nmax)->default_value("64"))
   ("m,mmax", "maximum cylindrical azimuthal order",
     cxxopts::value<int>(mmax)->default_value("6"))
   ("N,norder", "maximum cylindrical radial order",
     cxxopts::value<int>(norder)->default_value("18"))
   ("d,nodd", "number of vertically antisymmetric functions per M-order",
     cxxopts::value<int>(nodd)->default_value("-1"))
   ("a,ascale", "disk scale length",
     cxxopts::value<double>(acyl)->default_value("0.01"))
   ("h,hscale", "disk scale height",
     cxxopts::value<double>(hcyl)->default_value("0.001"))
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
  EmpCylSL::CMAPR       = 1;
  EmpCylSL::CMAPZ       = 1;
  EmpCylSL::DENS        = true;
  EmpCylSL::VFLAG       = 26;

  EmpCylSL test1(nmax, lmax, mmax, norder, acyl, hcyl, nodd, eof1);

  bool cache_ok = test1.read_cache();

  if (!cache_ok) {
    if (myid==0) {		// Diagnostic output . . .
      std::cerr << "Can not read explicitly specified EOF file <"
		<< eof1 << ">" << std::endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  EmpCylSL test2(nmax, lmax, mmax, norder, acyl, hcyl, nodd, eof2);

  cache_ok = test2.read_cache();

  if (!cache_ok) {
    if (myid==0) {		// Diagnostic output . . .
      std::cerr << "Can not read explicitly specified EOF file <"
		<< eof2 << ">" << std::endl;
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

