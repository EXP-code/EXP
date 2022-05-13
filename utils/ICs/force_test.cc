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

#include <DiskEval.H>

#include <fenv.h>

#include <Progress.H>		// Progress bar
#include <cxxopts.H>		// Option parsing

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
  
  cxxopts::Options options(av[0], "Force accuracy testing for cylindrical basis");

  options.add_options()
   ("h,help", "Print this help message")
   ("L,logr", "Use log grid for DiskEval")
   ("dmodel", "Target model type (MN or exponential)",
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
    const char *mesg = "Force errors using DiskEval";
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

  // Open mass, position [3], acceleration [3] data
  //
  std::ifstream in(fdata);
  if (!in) {
    std::cout << "Error opening <" << fdata << ">" << std::endl;
    exit(-2);
  }
  int nbods = 0;
  in.read((char *)&nbods, sizeof(int));
    
  // Make sure that nout is even to prevent divide by zero
  //
  if (2*(nout/2) != nout) nout++;

  double dR = (rout - rinn)/(nout-1);
  double dZ = 2.0*zout/(nout-1);

  using Dvector = std::vector< std::vector<double> >;

  Dvector mass(nout), meanFR2(nout), meanFz2(nout);

  for (auto & v : mass)    v.resize(nout);
  for (auto & v : meanFR2) v.resize(nout);
  for (auto & v : meanFz2) v.resize(nout);

  std::cout << std::endl << "Begin: particle force eval   "
	    << std::endl << "-----------------------------"
	    << std::endl;
  
  progress::progress_display progress(nbods);

  while (true) {

    ++progress;

    if (in.good()) {
      float m, pos[3], acc[3];

      in.read((char *)&m, sizeof(float));
      in.read((char *)&pos[0], sizeof(float)*3);
      in.read((char *)&acc[0], sizeof(float)*3);
      if (not in.good()) break;

      double R  = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      double z  = pos[2];

      double fR = (pos[0]*acc[0] + pos[1]*acc[1])/(R+1.0e-18);
      double fz = acc[2];

      if (R>=rinn and R<rout and z>=-zout and z<zout) {
	int iR = (R - rinn)/dR;
	int iZ = (z + zout)/dZ;

	iR = std::min<int>(iR, nout-1);
	iZ = std::min<int>(iZ, nout-1);

	auto ret = test(R, z);

	double fR_0 = std::get<1>(ret) * dmass;
	double fz_0 = std::get<2>(ret) * dmass;

	double dif  = fR - fR_0;

	mass   [iR][iZ] += m;
	meanFR2[iR][iZ] += m * (fR - fR_0) * (fR - fR_0);
	meanFz2[iR][iZ] += m * (fz - fz_0) * (fz - fz_0);
      }
    } else {
      break;
    }
  }

  // Plot potential and force plane evaluation in gnuplot format
  //
  std::ofstream out(outfile);

  if (out) {
    
    // Headers
    out << "#"
	<< std::setw(15) << std::right << "R |"
	<< std::setw(16) << std::right << "z |"
	<< std::setw(16) << std::right << "D(f_R)/f_R |"
	<< std::setw(16) << std::right << "D(f_z)/f_z |"
	<< std::setw(16) << std::right << "D(f_R)/scl |"
	<< std::setw(16) << std::right << "D(f_z)/scl |"
	<< std::endl
	<< "#"
	<< std::setw(15) << std::right << "[1] |"
	<< std::setw(16) << std::right << "[2] |"
	<< std::setw(16) << std::right << "[3] |"
	<< std::setw(16) << std::right << "[4] |"
	<< std::setw(16) << std::right << "[5] |"
	<< std::setw(16) << std::right << "[6] |"
	<< std::endl
	<< "#" << std::setfill('-')
	<< std::setw(15) << std::right << "+"
	<< std::setw(16) << std::right << "+"
	<< std::setw(16) << std::right << "+"
	<< std::setw(16) << std::right << "+"
	<< std::setw(16) << std::right << "+"
	<< std::setw(16) << std::right << "+"
	<< std::endl << std::setfill(' ');

    std::cout << std::endl << "Begin: force bin evaluation  "
	      << std::endl << "-----------------------------"
	      << std::endl;
      
    progress::progress_display progress(nout*nout);
  

    auto ret = test(A, 0.0);
    double fR_0 = std::fabs(std::get<1>(ret)) * dmass;

    for (int i=0; i<nout; i++) {
      double R = rinn + dR*(0.5+i);

      auto ret = test(R, H);

      double fz_0 = std::fabs(std::get<2>(ret)) * dmass;

      for (int j=0; j<nout; j++) {
	double z = -zout + dZ*(0.5+j);

	double ms = mass[i][j] + 1.0e-18;
	double stdFR = std::sqrt(meanFR2[i][j]/ms);
	double stdFZ = std::sqrt(meanFz2[i][j]/ms);

	auto ret1 = test(R, z);

	out << std::setw(16) << R
	    << std::setw(16) << z
	    << std::setw(16) << stdFR/fabs(std::get<1>(ret1)*dmass)
	    << std::setw(16) << stdFZ/fabs(std::get<2>(ret1)*dmass)
	    << std::setw(16) << stdFR/fR_0
	    << std::setw(16) << stdFZ/fz_0
	    << std::endl;

	++progress;
      }
      out << std::endl;
    }
    
  } else {
    std::cout << "Error opening output file" << std::endl;
  }

  return 0;
}

