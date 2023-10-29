/*
  A simple cube IC generator for normally distributed velocities,
  optionally anisotropic
*/
                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include <cxxopts.H>

int 
main(int ac, char **av)
{
  //=====================
  // Begin option parsing
  //=====================

  int          N;		// Number of particles
  double       M;		// Total mass
  std::string  bodyfile;	// Output file
  unsigned     seed;		// Will be inialized by /dev/random if
				// not set on the command line

  // Default values for the velocity dispersion and bulk velocity
  //
  std::vector<double> disp(3, 1.0), bulk(3, 0.0);
  
  cxxopts::Options options(av[0], "Cube IC generator with uniform spatial and normal velocities.  You may specify a bulk velocity and anisotropic dispersion vector.");

  options.add_options()
    ("h,help", "Print this help message")
    ("N,number", "Number of particles to generate",
     cxxopts::value<int>(N)->default_value("100000"))
    ("M,mass", "Total mass of the cube",
     cxxopts::value<double>(M)->default_value("1.0"))
    ("d,disp", "Velocity dispersion triple",
     cxxopts::value<std::vector<double>>(disp))
    ("v,bulk", "Bulk velocity triple",
     cxxopts::value<std::vector<double>>(bulk))
    ("s,seed", "Random number seed. Default: use /dev/random",
     cxxopts::value<unsigned>(seed))
    ("o,file", "Output body file",
     cxxopts::value<std::string>(bodyfile)->default_value("cube.bods"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    std::cout << options.help() << std::endl << std::endl;
    return 1;
  }

  // Set from /dev/random if not specified
  if (vm.count("seed")==0) {
    seed = std::random_device{}();
  }

  // Make sure N>0
  if (N<=0) {
    std::cerr << av[0] << ": you must requiest at least one body"
	      << std::endl;
  }

  // Open the output file
  //
  std::ofstream out(bodyfile);

  if (out) {
    // Header
    //
    out << std::setw(10) << N << std::setw(10) << 0 << std::setw(10) << 0
	<< std::endl << std::setprecision(10);

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    std::normal_distribution<> normal(0.0, 1.0);

    double mass = M/N;

    for (int n=0; n<N; n++) {
      out << std::setw(18) << mass;
      for (int k=0; k<3; k++) out << std::setw(18) << uniform(gen);
      for (int k=0; k<3; k++) out << std::setw(18) << normal(gen)*disp[k] + bulk[k];
      out << std::endl;
    }

  } else {
    std::cout << "Error opening output body file <"
	      << bodyfile << ">" << std::endl;
  }
  

  return 0;
}
