/*
  Use DiskEval to compute force errors from a mass, position, acceleration grid
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

#include "DiskEval.H"

#include <fenv.h>

#include "Progress.H"		// Progress bar
#include "cxxopts.H"		// Command line parsing

// Globals for exp libraries
//
#include "libvars.H"

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
  int          mmax;
  int          nump;
  double       dmass;
  double       A;
  double       H;
  double       rmin;
  double       rmax;
  double       rinn;
  double       rout;
  double       zout;
  string       dmodel;
  string       fdata;
  string       outfile;
  
  std::ostringstream sout;
  sout << "Use DiskEval to compute force errors from a mass, position, acceleration grid." << std::endl << "Evaluated at specific input locations." << std::endl;

  cxxopts::Options options(av[0], sout.str());

  options.add_options()
    ("h,help", "Print this help message")
    ("L,logr", "Use log grid for DiskEval")
    ("dmodel", "Target model type (Ferrers, MN, or exponential)",
     cxxopts::value<std::string>(dmodel)->default_value("exponential"))
    ("force", "Force data from N-body evluation",
     cxxopts::value<std::string>(fdata)->default_value("force.data"))
    ("out", "Output force test grid data",
     cxxopts::value<std::string>(outfile)->default_value("testforce.dat"))
    ("dmass", "Total disk mass",
     cxxopts::value<double>(dmass)->default_value("0.025"))
    ("nint", "Number of Gauss-Legendre knots for theta integration",
     cxxopts::value<int>(nint)->default_value("40"))
    ("numr", "Size of radial grid",
     cxxopts::value<int>(numr)->default_value("1000"))
    ("nump", "Size of azimuthal grid",
     cxxopts::value<int>(nump)->default_value("1"))
    ("lmax", "Number of harmonics for Spherical SL for halo/spheroid",
     cxxopts::value<int>(lmax)->default_value("32"))
    ("mmax", "Maximum m order",
     cxxopts::value<int>(mmax)->default_value("0"))
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
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
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

  
  bool logr = false;
  if (vm.count("logr")) logr = true;

  // Sanity check for inner logarithmic radius
  //
  if (rinn<=0.0) rinn = 1.0e-4;

  // The model instance (you can add others in DiskModels.H)
  //
  EmpCylSL::AxiDiskPtr model;

  if (dmodel.compare("FEED")==0) {// Ferrers + Evacuated Exponential Disk. most options are hardwired here relative to the disc scale length currently
    std::cout << "Using model FEED" << endl;
    model = std::make_shared<FEED>(1.4*A,(0.4/1.4)*A, (0.2/1.4)*A, A, H, 0.2);
  } else if (dmodel.compare("Ferrers")==0) {// Ferrers
    std::cout << "Using model Ferrers" << endl;
    model = std::make_shared<Ferrers>(A, H, A/10.);
  } else if (dmodel.compare("MN")==0) {// Miyamoto-Nagai
    std::cout << "Using Miyamoto-Nagai disk model" << endl;
    model = std::make_shared<MNdisk>(A, H);
  } else if (dmodel.compare("DoubleExponential")==0) {// Double Exponential. Most options are hardwired right now.
    std::cout << "Using Double Exponential disk model" << endl;
    model = std::make_shared<DoubleExponential>(A, H, A, H/3, 0.5);
  } else {			// Default to exponential
    std::cout << "Using standard exponential model" << endl;
    model = std::make_shared<Exponential>(A, H);
  }

  DiskEval test(model, rmin, rmax, A, lmax, numr, nint, true, mmax, nump);

  // Open mass, position [3], acceleration [3] data
  //
  std::ifstream in(fdata);
  if (!in) {
    std::cout << "Error opening <" << fdata << ">" << std::endl;
    exit(-2);
  }

  // this is BAD hardwire right now, need to read from fdata
  int nbods = 1000000;
  //in.read((char *)&nbods, sizeof(int));
    

  std::cout << std::endl << "Begin: particle force eval   "
	    << std::endl << "-----------------------------"
	    << std::endl;

  
  
  progress::progress_display progress(nbods);

  std::ofstream out(outfile);

  
  while (true) {

    ++progress;

    if (in.good()) {
      float m, x,y,z,ax,ay,az;

      in.read((char *)&m, sizeof(float));
      in.read((char *)&x, sizeof(float));
      in.read((char *)&y, sizeof(float));
      in.read((char *)&z, sizeof(float));
      in.read((char *)&ax, sizeof(float));
      in.read((char *)&ay, sizeof(float));
      in.read((char *)&az, sizeof(float));
      if (not in.good()) break;

      double R  = std::sqrt(x*x + y*y);
      double phi = atan2(y,x);

      double fR = (x*ax + y+ay)/(R+1.e-18);
      float fz = az;

      auto ret = test(R, z, phi);

      float potval = std::get<0>(ret) * dmass;
      float fR_0   = std::get<1>(ret) * dmass;
      float fz_0   = std::get<2>(ret) * dmass;
      float fphi_0 = std::get<3>(ret) * dmass;


      //cout << R << endl;
      //cout << fR_0 << endl;

      float f;
      out.write((const char *)&(f=x), sizeof(float));
      out.write((const char *)&(f=y), sizeof(float));
      out.write((const char *)&(f=z), sizeof(float));
      out.write((const char *)&(f=fR_0), sizeof(float));
      out.write((const char *)&(f=fz_0), sizeof(float));
      out.write((const char *)&(f=fphi_0), sizeof(float));
      out.write((const char *)&(f=potval), sizeof(float));


    }
  }


  return 0;
}


