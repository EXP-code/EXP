//  Generates a cylindrical basis cache file using expui

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <config_exp.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

// EXP classes
//
#include <libvars.H>
#include <BasisFactory.H>
#include <cxxopts.H>		// Command-line parsing
#include <EXPini.H>		// Ini-style config

int
main(int ac, char **av)
{
  // Parameters
  //
  int nthrds = 1;
  std::string config;

  // Default/example YAML file for disk
  //
  std::string disk_config = 
    "id           : cylinder\n"
    "parameters   :\n"
    "  acyl       : 1.0\n"
    "  hcyl       : 0.1\n"
    "  lmaxfid    : 48\n"
    "  nmaxfid    : 48\n"
    "  mmax       : 10\n"
    "  nmax       : 32\n"
    "  ncylodd    : 6\n"
    "  ncylnx     : 256\n"
    "  ncylny     : 128\n"
    "  ncylr      : 2000\n"
    "  rnum       : 200\n"
    "  pnum       : 1\n"
    "  tnum       : 80\n"
    "  rcylmin    : 0.001\n"
    "  rcylmax    : 20\n"
    "  ashift     : 0\n"
    "  logr       : true\n"
    "  expcond    : true\n"
    "  deproject  : true\n"
    "  cachename  : eof.cache.file_new\n"
    "  vflag      : 5\n"
    ;

  // Inialize MPI stuff
  //
  local_init_mpi(ac, av);

  const std::string mesg("Generates an EmpCylSL cache from a YAML configuration file\n");

  cxxopts::Options options(av[0], mesg);

  options.add_options()
    ("h,help", "Print this help message")
    ("T,template", "Write template YAML config file")
    ("c,config", "Read a configuration file", cxxopts::value<std::string>(config))
    ("t,threads", "Number of OpenMP threads", cxxopts::value<int>(nthrds))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Write YAML template config file and exit
  //
  if (vm.count("template")) {
    if (myid==0) {
      std::ofstream out("template.yaml");
      if (out) {
	out << disk_config;
	std::cout << "Wrote <template.yaml>.  Please adjust to your preferences"
		  << std::endl;
      } else {
	std::cerr << "Error opening <template.yaml> for output" << std::endl;
      }
    }
    MPI_Finalize();
    return 0;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << options.help() << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the YAML config file
  //
  if (vm.count("config")) {
    try {
      auto name = vm["config"].as<std::string>();
      auto size = std::filesystem::file_size(name);
      std::string content(size, '\0');
      std::ifstream in(name);
      in.read(&content[0], size);
      disk_config = content;
    } catch (std::runtime_error& e) {
      if (myid==0) std::cout << "Error reading configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }

  // Explicit OpenMP control; otherwise, fall back to the environment
  // setting OMP_NUM_THREADS.
  //
#ifdef HAVE_OMP_H
  if (vm.count("threads")) omp_set_num_threads(nthrds);
#pragma omp parallel
  if (myid==0) {
    int numthrd = omp_get_num_threads();
    std::cout << "Number of threads=" << numthrd << std::endl;
  }
#endif

   // Begin calculation and start stopwatch
   //
   auto t_start = MPI_Wtime();

   // Construct the basis instance
   //
   auto disk_basis = BasisClasses::Basis::factory_string(disk_config);

   MPI_Barrier(MPI_COMM_WORLD);

   // DONE
   //
   if (myid==0)
     std::cout << std::string(60, '=')  << std::endl
	       << "Calculation finished in " << std::setprecision(3)
	       << MPI_Wtime() - t_start << " seconds" << std::endl
	       << std::string(60, '=') << std::endl;

  // Shutdown MPI
  //
  MPI_Finalize();

  return 0;
}
