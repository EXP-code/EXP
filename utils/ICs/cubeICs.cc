/*
  A simple cube IC generator for normally distributed velocities,
  optionally anisotropic.  Also includes a perturbation wave vector
  test, generated by acceptance-rejection sampling.
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
  std::vector<double> disp, bulk;

  // Perturbation wave vector and amplitude
  //
  std::vector<int> pert;
  double pamp = 0.0;
  
  cxxopts::Options options(av[0], "Cube IC generator with uniform spatial and normal velocities.  You may specify a bulk velocity and anisotropic dispersion vector.");

  options.add_options()
    ("h,help",    "Print this help message")
    ("z,zerovel", "Zero the mean velocity")
    ("N,number",  "Number of particles to generate",
     cxxopts::value<int>(N)->default_value("100000"))
    ("M,mass",    "Total mass of the cube",
     cxxopts::value<double>(M)->default_value("1.0"))
    ("d,disp",    "Velocity dispersion triple",
     cxxopts::value<std::vector<double>>(disp))
    ("v,bulk",    "Bulk velocity triple",
     cxxopts::value<std::vector<double>>(bulk))
    ("s,seed",    "Random number seed. Default: use /dev/random",
     cxxopts::value<unsigned>(seed))
    ("o,file",    "Output body file",
     cxxopts::value<std::string>(bodyfile)->default_value("cube.bods"))
    ("w,wave",    "Perturbation wave vector",
     cxxopts::value<std::vector<int>>(pert))
    ("p,pamp",    "Perturbation amplitude",
     cxxopts::value<double>(pamp))
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

  // Check for consistent pertubation amplitude
  bool pwave = false;
  if (vm.count("pamp")) {
    if (pert.size()!=3) {
      std::cout << "You must specify a three vector wave vector"
		<< ". Current size=" << pert.size()
		<< std::endl;
      return 2;
    }
    int ntot = 0;
    for (int k=0; k<3; k++) ntot += pert[k]*pert[k];
    if (ntot==0) {
      std::cout << "You must specify a non-zero pertubation wave number"
		<< std::endl;
      return 2;
    }
    if (fabs(pamp)>=1.0) {
      std::cout << "Your pertubation amplitude must be in (-1, 1)"
		<< std::endl;
      return 2;
    }
    pwave = true;
  }

  if (disp.size()<3) {
    for (auto j=disp.size(); j<3; j++) disp.push_back(1.0);
  }

  if (bulk.size()<3) {
    for (auto j=bulk.size(); j<3; j++) bulk.push_back(0.0);
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

    // Save the position and velocity vectors
    std::vector<std::array<double, 3>> pos, vel;

    // Temporary position for iteration
    std::array<double, 3> x;

    // Generation loop
    for (int n=0; n<N; n++) {

      // Select a wave perturbation?
      if (pwave) {

	// Acceptance-rejection loop
	while (true) {
	  double nx = 0.0;
	  for (int k=0; k<3; k++) {
	    x[k] = uniform(gen);
	    nx  += x[k] * pert[k];
	  }
	  double frac = (1.0 + pamp*cos(2.0*M_PI*nx))/(1.0 + fabs(pamp));
	  if (frac > uniform(gen)) break; // Accept?
	}
	// Make the position vector
	pos.push_back(x);
      }
      else {
	// Make the position vector
	pos.push_back({uniform(gen), uniform(gen), uniform(gen)});
      }

      // Make the velocity vector
      for (int k=0; k<3; k++) x[k] = normal(gen)*disp[k] + bulk[k];
      vel.push_back(x);
    }

    // Zero the mean velocity
    //
    if (vm.count("zerovel")) {
      x = {0, 0, 0};
      double w = 1.0/N;
      for (int i=0; i<N; i++) {
	for (int k=0; k<3; k++) x[k] += w*vel[i][k];
      }
      for (int i=0; i<N; i++) {
	for (int k=0; k<3; k++) vel[i][k] -= x[k];
      }
    }
      
    double mass = M/N;

    for (int i=0; i<N; i++)  {
      out << std::setw(18) << mass;
      for (int k=0; k<3; k++) out << std::setw(18) << pos[i][k];
      for (int k=0; k<3; k++) out << std::setw(18) << vel[i][k];
      out << std::endl;
    }

  } else {
    std::cout << "Error opening output body file <"
	      << bodyfile << ">" << std::endl;
  }
  

  return 0;
}
