/*
  Add a ring of particles to an existing realization
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

#include <Eigen/Eigen>

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
#include <hernquist_model.H>
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

  std::string  halofile;
  int          DIVERGE, N;
  double       DIVERGE_RFAC, mass, radius, width, vdisp;
  std::string  outfile, config;

  const std::string mesg("Make a model from the combination of two spherical models.  E.g. a DM halo and a bulge.\n");

  cxxopts::Options options(av[0], mesg);

  options.add_options()
    ("h,help", "Print this help message")
    ("s,spline", "Set interplation in SphericalModelTable to 'spline' (default: linear)")
    ("T,template", "Write template options file with current and all default values")
    ("c,config", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("N,number", "Number of particles",
     cxxopts::value<int>(N)->default_value("100000"))
    ("DIVERGE", "Cusp power-law density extrapolation (0 means off)",
     cxxopts::value<int>(DIVERGE)->default_value("0"))
    ("DIVERGE_RFAC", "Exponent for inner cusp extrapolation",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("1.0"))
    ("i,halofile", "Halo spherical model table",
     cxxopts::value<string>(halofile)->default_value("SLGridSph.one"))
    ("o,outfile", "Output particle file",
     cxxopts::value<string>(outfile)->default_value("addring.bods"))
    ("M,mass", "Mass of ring",
     cxxopts::value<double>(mass)->default_value("0.1"))
    ("R,radius", "Radius of ring",
     cxxopts::value<double>(radius)->default_value("1.0"))
    ("W,width", "Width of ring",
     cxxopts::value<double>(width)->default_value("0.1"))
    ("V,disp", "Velocity dispersion",
     cxxopts::value<double>(vdisp)->default_value("0.1"))
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

  // Open output file
  //
  std::ofstream out(outfile);
  if (not out) throw std::runtime_error("Cannot open output file " + outfile);

  // Model
  //
  SphericalModelTable model(halofile, DIVERGE, DIVERGE_RFAC);


  // Random setup
  //
  std::random_device rd{};
  std::mt19937 gen{rd()};
 
  // Unit normal distribution
  //
  std::normal_distribution<> d{0.0, 1.0};
  std::uniform_real_distribution<> u{0.0, 1.0};
  
 
  double pmass = mass/N;
  Eigen::Vector3d pos, vel;
  
  for (int i=0; i<N; i++) {

    double R  = radius + 0.5*width*d(gen);
    double v0 = sqrt(model.get_mass(radius)/radius);
    double phi = 2.0*M_PI*u(gen);

    pos(0) = R*cos(phi);
    pos(1) = R*sin(phi);
    pos(2) = 0.5*width*d(gen);
    
    vel(0) = vdisp*d(gen) - v0*sin(phi);
    vel(1) = vdisp*d(gen) + v0*cos(phi);
    vel(2) = vdisp*d(gen);
    
    out << std::setw(18) << pmass;
    for (int k=0; k<3; k++) out << std::setw(18) << pos(k);
    for (int k=0; k<3; k++) out << std::setw(18) << vel(k);
    out << std::endl;
  }

  return 0;
}
