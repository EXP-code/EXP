/*
  Add together two spherical models
*/

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <random>
#include <vector>
#include <cmath>

#include <fenv.h>

#include <config_exp.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

// EXP classes
//
#include <numerical.H>
#include <gaussQ.H>
#include <isothermal.H>
#include <hernquist.H>
#include <model3d.H>
#include <biorth.H>
#include <SphericalSL.H>
#include <interp.H>
#include <AddSpheres.H>
#include <libvars.H>		// Library globals
#include <cxxopts.H>		// Command-line parsing
#include <EXPini.H>		// Ini-style config

                                // Local headers
#include <SphericalSL.H>
#include <localmpi.H>

int 
main(int ac, char **av)
{
  //====================
  // Inialize MPI stuff
  //====================
  
  local_init_mpi(ac, av);

  //====================
  // Begin opt parsing
  //====================

  std::string  halofile1, halofile2;
  int          DIVERGE, DIVERGE2;
  double       DIVERGE_RFAC, DIVERGE_RFAC2, bmass;
  std::string  outfile, config;
  
  const std::string mesg("Make a model from the combination of two spherical models.  E.g. a DM halo and a bulge.\n");

  cxxopts::Options options(av[0], mesg);

  options.add_options()
    ("h,help", "Print this help message")
    ("s,spline", "Set interplation in SphericalModelTable to 'spline' (default: linear)")
    ("T,template", "Write template options file with current and all default values")
    ("c,config", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("DIVERGE", "Cusp power-law density extrapolation (0 means off)",
     cxxopts::value<int>(DIVERGE)->default_value("0"))
    ("DIVERGE_RFAC", "Exponent for inner cusp extrapolation",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("1.0"))
    ("DIVERGE2", "Cusp power-law extrapolation, second model (0 means off)",
     cxxopts::value<int>(DIVERGE2)->default_value("0"))
    ("DIVERGE_RFAC2", "Exponent for inner cusp extrapolation, second model",
     cxxopts::value<double>(DIVERGE_RFAC2)->default_value("1.0"))
    ("halofile1", "Halo spherical model table",
     cxxopts::value<string>(halofile1)->default_value("SLGridSph.one"))
    ("halofile2", "Second spherical model file",
     cxxopts::value<string>(halofile2)->default_value("SLGridSph.two"))
    ("newfile", "Output model name",
     cxxopts::value<string>(outfile)->default_value("new.model"))
    ("bmass", "Mass factor for second model relative to first",
     cxxopts::value<double>(bmass)->default_value("1.0"))
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
    // Write template file
    //
    if (myid==0) SaveConfig(vm, options, "template.yaml");

    MPI_Finalize();
    return 0;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << options.help() << std::endl << std::endl
		<< "Examples: " << std::endl
		<< "\t" << "Use parameters read from a YAML config file"  << std::endl
		<< "\t" << av[0] << " --config=addsphmod.config"  << std::endl << std::endl
		<< "\t" << "Generate a template YAML config file from current defaults"  << std::endl
		<< "\t" << av[0] << " --template=template.yaml" << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the YAML config file
  //
  if (vm.count("config")) {
    try {
      vm = LoadConfig(options, config);
    } catch (cxxopts::OptionException& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }
  
  if (vm.count("spline")) {
    SphericalModelTable::linear = 0;
  } else {
    SphericalModelTable::linear = 1;
  }

  auto mod1 = std::make_shared<SphericalModelTable>(halofile1, DIVERGE, DIVERGE_RFAC);

  auto mod2 = std::make_shared<SphericalModelTable>(halofile2, DIVERGE2, DIVERGE_RFAC2);


  AddSpheres combo(mod1, mod2, bmass);
  
  combo.get_model()->print_model(outfile);

  return 0;
}
