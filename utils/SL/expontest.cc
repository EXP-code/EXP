#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <Eigen/Eigen>

#include <PotRZ.H>		// Hankel computation for potential
#include <gaussQ.H>		// Gauss-Legendre quadrature
#include <cxxopts.H>

int main(int argc, char** argv)
{
  int M=0, N, nout;
  double A, rmax, rout;
  std::string filename;

  // Parse command line
  //
  cxxopts::Options options
    (argv[0],
     "Check PotRZ and QDHT using the exact exponential disk solution\n");
  
  options.add_options()
    ("h,help",     "Print this help message")
    ("N,nsize",    "Default radial grid size",
     cxxopts::value<int>(N)->default_value("256"))
    ("A,length",   "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("rmax",       "Outer radius for transform",
     cxxopts::value<double>(rmax)->default_value("10.0"))
    ("rout",       "Outer radius for evaluation",
     cxxopts::value<double>(rout)->default_value("10.0"))
    ("nout",       "number of points in the output grid per side",
     cxxopts::value<int>(nout)->default_value("40"))
    ("o,filename", "Output filename",
     cxxopts::value<std::string>(filename)->default_value("test.potrz"))
    ;

  
  //===================
  // Parse options
  //===================

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return 2;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    std::cout << options.help() << std::endl << std::endl;
    return 1;
  }

  // Open output file
  //
  std::ofstream out(filename);
  if (not out) throw std::runtime_error("Could not open output file: " +
					filename);

  // Define some representative limits
  //
  double Rmax = rout;

  // Grid spacing
  //
  double dR = Rmax/(nout - 1);

  // Potential instance with radially sensitive convergence parameters
  //
  PotRZ pot(rmax, N, M);

  // Set the functor using a lambda
  //
  auto dens = [A](double R) {
    return -exp(-R/A);
  };
      
  auto potl = [A](double R) {
    double x = 0.5*R/A;
    return M_PI*R*(std::cyl_bessel_i(1, x) * std::cyl_bessel_k(0, x) -
		   std::cyl_bessel_i(0, x) * std::cyl_bessel_k(1, x) );
  };

  // Write the results
  //
  for (int i=0; i<nout; i++) {
    double R = dR*i;
    out << std::setw(16) << R
	<< std::setw(16) << potl(R)
	<< std::setw(16) << pot(R, 0.0, dens)
	<< std::endl;
  }

  return 0;
}
