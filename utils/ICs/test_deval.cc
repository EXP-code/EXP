/*
  Check DiskEval implementation
*/
                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <DiskEval.H>
#include <cxxopts.H>

#include <fenv.h>

// Globals for exp libraries
//
#include <libvars.H>

int 
main(int ac, char **av)
{
  //====================
  // Begin opt parsing
  //====================

  int          lmax;
  int          numr;
  int          nint;
  int          nout;
  double       A;
  double       H;
  double       rmin;
  double       rmax;
  double       rinn;
  double       rout;
  double       zout;
  string       dmodel;
  
  cxxopts::Options options(av[0], "Check DiskEval implementation");

  options.add_options()
   ("h,help", "Print this help message")
   ("L,logr", "Use log grid for DiskEval")
   ("dmodel", "Target model type (MN or exponential)",
     cxxopts::value<std::string>(dmodel)->default_value("exponential"))
   ("nint", "Number of Gauss-Legendre knots for theta integration",
     cxxopts::value<int>(nint)->default_value("40"))
   ("numr", "Size of radial grid",
     cxxopts::value<int>(numr)->default_value("1000"))
   ("lmax", "Number of harmonics for Spherical SL for halo/spheroid",
     cxxopts::value<int>(lmax)->default_value("32"))
   ("rmin", "Minimum radius for grid",
     cxxopts::value<double>(rmin)->default_value("1.0e-4"))
   ("rmax", "Maximum radius for grid",
     cxxopts::value<double>(rmax)->default_value("1.0"))
   ("A", "Radial scale length for disk basis construction",
     cxxopts::value<double>(A)->default_value("0.01"))
   ("H", "Vertical scale length for disk basis construction",
     cxxopts::value<double>(H)->default_value("0.001"))
   ("r", "Minimum cylindrical radius for test output",
     cxxopts::value<double>(rinn)->default_value("0.0001"))
   ("R", "Maximum cylindrical radius for test output",
     cxxopts::value<double>(rout)->default_value("0.3"))
   ("Z", "Maximum height for testoutput",
     cxxopts::value<double>(zout)->default_value("0.1"))
   ("N", "Number of grid points for test plane",
     cxxopts::value<int>(nout)->default_value("60"))
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

  
  bool logr = false;
  if (vm.count("logr")) logr = true;

  // Sanity check for inner logarithmic radius
  //
  if (rinn<=0.0) rinn = 1.0e-4;

  // The model instance (you can add others in DiskModels.H)
  //
  EmpCylSL::AxiDiskPtr model;
      
  if (dmodel.compare("MN")==0) // Miyamoto-Nagai
    model = std::make_shared<MNdisk>(A, H);
  else			// Default to exponential
    model = std::make_shared<Exponential>(A, H);
      
  DiskEval test(model, rmin, rmax, A, lmax, numr, nint, true);

  // Plot potential and force plane evaluation in gnuplot format
  //
  std::ofstream out("testdeval.eval");

  if (out) {
				// Logarithmic output grid in radius
    double dR = (log(rout) - log(rinn))/(nout - 1);
    double dz = 2.0*zout/(nout - 1);

    for (int i=0; i<nout; i++) {
      double R = rinn*exp(dR*i);

      for (int j=0; j<nout; j++) {
	double z = -zout + dz*j;

	auto ret = test(R, z);

	// Miyamoto-Nagai test
	//
	double pot_t = 0.0, FR_t = 0.0, Fz_t = 0.0;
	if (dmodel.compare("MN")==0) {
	  double zb = sqrt(z*z + H*H);
	  double s  = A + zb;
	  pot_t = -1.0/sqrt(R*R + s*s);
	  FR_t  = -R*pow(R*R + s*s, -1.5);
	  Fz_t  = -z/zb*s*pow(R*R + s*s, -1.5);
	}

	out << std::setw(18) << R
	    << std::setw(18) << z
	    << std::setw(18) << std::get<0>(ret)
	    << std::setw(18) << std::get<1>(ret)
	    << std::setw(18) << std::get<2>(ret)
	    << std::setw(18) << pot_t
	    << std::setw(18) << FR_t
	    << std::setw(18) << Fz_t
	    << std::endl;
      }
      out << std::endl;
    }

  } else {
    std::cout << "Error opening output file" << std::endl;
  }
  
  out.close();
  out.open("testdeval.midplane");

  if (out) {
    double dR = (log(rout) - log(rinn))/(nout - 1);
    
    for (int i=0; i<nout; i++) {
      double R = rinn*exp(dR*i);

      auto ret = test(R, 0);

      out << std::setw(18) << R
	  << std::setw(18) << std::get<0>(ret)
	  << std::setw(18) << std::get<1>(ret)
	  << std::setw(18) << std::get<2>(ret)
	  << std::endl;
    }
    
  } else {
    std::cout << "Error opening output file" << std::endl;
  }
  

  return 0;
}

