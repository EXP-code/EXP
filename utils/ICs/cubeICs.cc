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
#include <array>

#include <cxxopts.H>

int 
main(int ac, char **av)
{
  //====================
  // Begin opt parsing
  //====================

  int          N;		// Number of particles
  double       M;		// Total mass
  unsigned     seed;
  std::array<double, 3> disp;
  std::string  bodyfile;
  
  cxxopts::Options options(av[0], "Cube IC generator");

  options.add_options()
    ("h,help", "Print this help message")
    ("N,number", "Number of particles to generate",
     cxxopts::value<int>(N)->default_value("100000"))
    ("M,mass", "Total mass of the cube",
     cxxopts::value<double>(M)->default_value("1.0"))
    ("x,dispX", "Velocity disperion in the X direction",
     cxxopts::value<double>(disp[0])->default_value("1.0"))
    ("y,dispY", "Velocity disperion in the Y direction",
     cxxopts::value<double>(disp[1])->default_value("1.0"))
    ("z,dispZ", "Velocity disperion in the Z direction",
     cxxopts::value<double>(disp[2])->default_value("1.0"))
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

  if (vm.count("seed")==0) {
    seed = std::random_device{}();
  }

  if (N<=0) {
    std::cerr << argv[0] << ": you must requiest at least one body"
	      << std::endl;
  }

  // Open the oputput file
  //
  std::ofstream out(bodyfile);

  if (out) {
    // Header
    //
    out << std::setw(10) << N << std::setw(10) << 0 << std::setw(10)
	<< std::endl;
    
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    std::normal_distribution<> normal(0.0, 1.0);

    double mass = Mass/N;

    for (int n=0; n<N; n++) {
      out << std::setw(18) << mass;
      for (int k=0; k<3; k++) out << std::setw(18) << uniform(gen);
      for (int k=0; k<3; k++) out << std::setw(18) << normal(gen)*disp[k];
      out << std::endl;
    }

  } else {
    std::cout << "Error opening output body file <"
	      << bodyfile << ">" << std::endl;
  }
  

  return 0;
}

