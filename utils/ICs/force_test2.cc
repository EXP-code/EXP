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

// Boost stuff
//
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>	// Progress bar

// Globals for exp libraries
//
int myid = 0;
char threading_on  = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
string outdir, runtag;


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
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("logr,L",                                                                          "Use log grid for DiskEval")
    ("dmodel",          po::value<std::string>(&dmodel)->default_value("exponential"),  "Target model type (Ferrers, MN or exponential)")
    ("force",           po::value<std::string>(&fdata)->default_value("force.data"),  "Force data from N-body evluation")
    ("out",             po::value<std::string>(&outfile)->default_value("testforce.dat"),  "Output force test grid data")
    ("dmass",           po::value<double>(&dmass)->default_value(0.025),  "Total disk mass")
    ("nint",            po::value<int>(&nint)->default_value(40),                       "Number of Gauss-Legendre knots for theta integration")
    ("numr",            po::value<int>(&numr)->default_value(1000),                     "Size of radial grid")
    ("nump",            po::value<int>(&nump)->default_value(1), "Size of azimuthal grid")
    ("lmax",            po::value<int>(&lmax)->default_value(32),                       "Number of harmonics for Spherical SL for halo/spheroid")
    ("mmax",            po::value<int>(&mmax)->default_value(0), "Maximum m order")
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
    const char *mesg = "Force errors using DiskEval";
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

  if (dmodel.compare("FEED")==0) // Ferrers + Evacuated Exponential Disk. most options are hardwired here relative to the disc scale length currently
    model = boost::make_shared<FEED>(1.4*A,(0.4/1.4)*A, (0.2/1.4)*A, A, H, 0.2);
  else if (dmodel.compare("Ferrers")==0) // Ferrers
    model = boost::make_shared<Ferrers>(A, H, A/10.);
  else if (dmodel.compare("MN")==0) // Miyamoto-Nagai
    model = boost::make_shared<MNdisk>(A, H);
  else			// Default to exponential
    model = boost::make_shared<Exponential>(A, H);
      
  DiskEval test(model, rmin, rmax, A, lmax, numr, nint, true, mmax, nump);

  // Open mass, position [3], acceleration [3] data
  //
  std::ifstream in(fdata);
  if (!in) {
    std::cout << "Error opening <" << fdata << ">" << std::endl;
    exit(-2);
  }
  int nbods = 1000000;
  //in.read((char *)&nbods, sizeof(int));
    
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

  
  
  boost::progress_display progress(nbods);

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


