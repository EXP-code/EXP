#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <EmpCyl2D.H>
#include <cxxopts.H>

#include "PotRZ.H"

int main(int argc, char** argv)
{
  bool logr = false, cmap = false, ortho = false;
  int numr, mmax, nmax, knots, M, N;
  double A, scale, rmin, rmax;
  std::string filename, type, biorth;

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
    ("vertical", "Compute the vertical grid")
    ("scale", "scaling from real coordinates to table",
     cxxopts::value<double>(scale)->default_value("1.0"))
    ("M,harmonic", "Aximuthal harmonic m=0,1,2,3,...",
     cxxopts::value<int>(M)->default_value("0"))
    ("N,norder", "Radial harmonic for rendering",
     cxxopts::value<int>(N)->default_value("1"))
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
    ("biorth", "Biorthogonal type (cb, bess)",
     cxxopts::value<std::string>(biorth)->default_value("cb"))
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
  EmpCyl2D emp(mmax, nmax, knots, numr, rmin, rmax, A, 1.0, cmap, logr,
	       type, biorth);

  // Sanity check
  //
  M = std::min<int>(M, mmax);

  if (vm.count("grid"))  emp.writeBasis(M, filename + ".grid");
  if (vm.count("trans")) emp.writeTrans(M, filename + ".trans");
  if (vm.count("ortho")) emp.orthoCheck(M, filename + ".ortho");

  emp.checkCoefs();

  if (vm.count("vertical")) {

    // Create the functor
    auto mass = [&emp, M, N](double x)
    {
      return x * emp.get_dens(x, M, N);
    };

    // Grid size
    //
    constexpr int num = 40;

    // Define some representative limits
    //
    double Rmin = 0.001*A;
    double Rmax = 3.00*A;
    double Zmax = 1.00*A;

    // Grid spacing
    //
    double dz = Zmax/(num - 1);
    double dr = (Rmax - Rmin)/(num - 1);
  
    Eigen::MatrixXd outP(num, num), outH(num, num);
    Eigen::MatrixXi outN(num, num);

    // Do the grid computation
    //
    for (int j=0; j<num; j++) {
	
      double r = Rmin + dr*j;
      
      // These quadrature parameters are all empirical based on the
      // exponential disk

      // For inverse convergence
      //
      double h1 = std::min<double>(0.003*r, 0.05);
      int    N1 = std::min<int>(1000, floor(std::max<double>(100, 0.5/h1)));

      // For Sk convergence
      //
      double h2 = 0.05;
      int    N2 = 60;
      
      // Potential instance with radially sensitive convergence parameters
      //
      PotRZ pot(h1, h2, N1, N2, M, mass);

      for (int i=0; i<num; i++) {
	outP(j, i) = pot(r, dz*i);
	outH(j, i) = h1;
	outN(j, i) = N1;
      }
    }

    // Output file for grid
    //
    std::ofstream out(filename + ".RZ");

    out << std::setw(8) << num << std::setw(8) << num << std::endl;

    // Print the grid in the expected order
    //
    for (int i=0; i<num; i++) {
      for (int j=0; j<num; j++) {
	out << std::setw(16) << Rmin + dr*j
	    << std::setw(16) << dz*i
	    << std::setw(16) << outP(j, i)
	    << std::setw(16) << outH(j, i)
	    << std::setw( 8) << outN(j, i)
	    << std::endl;
      }
    }
  }

  return 0;
}
