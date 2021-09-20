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

#include <fenv.h>

// Boost stuff
//
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// Globals for exp libraries
//
#include <global.H>

namespace po = boost::program_options;

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
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("logr,L",                                                                          "Use log grid for DiskEval")
    ("dmodel",          po::value<std::string>(&dmodel)->default_value("exponential"),  "Target model type (MN or exponential)")
    ("nint",            po::value<int>(&nint)->default_value(40),                       "Number of Gauss-Legendre knots for theta integration")
    ("numr",            po::value<int>(&numr)->default_value(1000),                     "Size of radial grid")
    ("lmax",            po::value<int>(&lmax)->default_value(32),                       "Number of harmonics for Spherical SL for halo/spheroid")
    ("rmin",            po::value<double>(&rmin)->default_value(1.0e-4),                 "Minimum radius for grid")
    ("rmax",            po::value<double>(&rmax)->default_value(1.0),                    "Maximum radius for grid")
    ("A",               po::value<double>(&A)->default_value(0.01),                      "Radial scale length for disk basis construction")
    ("H",               po::value<double>(&H)->default_value(0.001),                     "Vertical scale length for disk basis construction")
    ("r",               po::value<double>(&rinn)->default_value(0.0001),                 "Minimum cylindrical radius for test output")
    ("R",               po::value<double>(&rout)->default_value(0.3),                    "Maximum cylindrical radius for test output")
    ("Z",               po::value<double>(&zout)->default_value(0.1),                    "Maximum height for testoutput") 
    ("N",               po::value<int>(&nout)->default_value(60),                        "Number of grid points for test plane");
       
  po::variables_map vm;
  
  // Parse command line for control and critical parameters
  //
  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error on command line: "
	      << e.what() << std::endl;
    return -1;
  }
  
  // Print help message and exit
  //
  if (vm.count("help")) {
    const char *mesg = "Test DiskEval implementation";
    std::cout << mesg << std::endl
	      << desc << std::endl << std::endl;
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
    model = boost::make_shared<MNdisk>(A, H);
  else			// Default to exponential
    model = boost::make_shared<Exponential>(A, H);
      
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

