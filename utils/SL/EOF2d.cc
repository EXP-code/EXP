#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <EmpCyl2D.H>
#include <cxxopts.H>

int main(int argc, char** argv)
{
  bool logr = false, cmap = false, ortho = false;
  int numr, mmax, nmax, knots, M;
  double A, scale, rmin, rmax;
  std::string filename, type;

  // Parse command line
  //
  cxxopts::Options options(argv[0],
			   "Computes an EOF two-dimensional disk basis from the Clutton-Brock basis for\n"
			   "one of the Kuzmin, finite Mestel or Exponential disk targets.  [The Mestel\n"
			   "disk will work very poorly because the Clutton-Brock basis has infinite\n"
			   "support and looks nothing like the Mestel disk profile.] The new basis, the\n"
			   "orthgonogality matrix and the rotation matrices may be written to files.\n");
  options.add_options()
    ("h,help", "Print this help message")
    ("logr", "Plot output grid with logarithmic spacing")
    ("cmap", "Use mapped coordinates")
    ("ortho", "Compute EOF orthogonal matrix and write to a file")
    ("grid", "Print the new basis grid to a file")
    ("trans", "Print the rotation matrices to a file")
    ("scale", "scaling from real coordinates to table",
     cxxopts::value<double>(scale)->default_value("1.0"))
    ("M,harmonic", "Aximuthal harmonic m=0,1,2,3,...",
     cxxopts::value<int>(M)->default_value("0"))
    ("A,length", "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("mmax", "maximum number of angular harmonics in the expansion",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("nmax", "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
    ("numr", "radial knots for the SL grid",
     cxxopts::value<int>(numr)->default_value("4000"))
    ("r,rmin", "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
    ("R,rmax", "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("20.0"))
    ("knots", "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("40"))
    ("type", "Target model type (kuzmin, mestel, expon)",
     cxxopts::value<std::string>(type)->default_value("expon"))
    ("o,filename", "Output filename",
     cxxopts::value<std::string>(filename)->default_value("testeof"))
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

  // Log spacing?
  //
  if (vm.count("logr")) logr = true;

  // Mapped coordinates?
  //
  if (vm.count("cmap")) cmap = true;

  // Make the class instance
  //
  EmpCyl2D emp(mmax, nmax, knots, numr, rmin, rmax, A, 1.0, true, logr, type);

  // Sanity check
  M = std::min<int>(M, mmax);

  if (vm.count("grid"))  emp.writeBasis(M, filename + ".grid");
  if (vm.count("trans")) emp.writeTrans(M, filename + ".trans");
  if (vm.count("ortho")) emp.orthoCheck(M, filename + ".ortho");

  return 0;
}
