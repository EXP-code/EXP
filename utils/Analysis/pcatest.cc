/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Dump pca basis
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 03/10/18
 *
 ***************************************************************************/

				// C++/STL headers
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <memory>

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// EXP classes
#include <global.H>
#include <numerical.H>
#include <interp.H>
#include <EmpCylSL.H>
#include <cxxopts.H>

#include <localmpi.H>
#include <foarray.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

const std::string overview = "Compute and print PCA basis for rendering";

// Globals
//

int
main(int argc, char **argv)
{
  double RMAX, ZMAX;
  int    OUTR, OUTZ, lmax, mmax, nmax, norder, numx, numy, nodd;
  std::string CACHEFILE, OUTTAG;
  double rcylmin, rcylmax, rscale, vscale;
  bool   DENS, LOGSC;

  //
  // Parse Command line
  //
  
  cxxopts::Options options("pcatest", "Unpack and summarize the pca basis");

  options.add_options()
    ("h,help", "produce this help message")
    ("v,verbose", "verbose output")
    ("R,RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("0.1"))
    ("Z,ZMAX", "maximum height for output",
     cxxopts::value<double>(ZMAX)->default_value("0.01"))
    ("rcylmin", "minimum radius for cylindrical basis table",
     cxxopts::value<double>(rcylmin)->default_value("0.001"))
    ("rcylmax", "maximum radius for cylindrical basis table",
     cxxopts::value<double>(rcylmax)->default_value("20.0"))
    ("rscale", "radial scale length for basis expansion",
     cxxopts::value<double>(rscale)->default_value("0.01"))
    ("vscale", "vertical scale length for basis expansion",
     cxxopts::value<double>(vscale)->default_value("0.001")) 
    ("lmax", "maximum harmonic order for spherical expansion",
     cxxopts::value<int>(lmax)->default_value("36"))
    ("nmax", "maximum harmonic order for spherical expansion",
     cxxopts::value<int>(nmax)->default_value("8"))
    ("mmax", "maximum azimuthal harmonic order for cylindrical expansion",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("norder", "maximum radial order for each harmonic subspace",
     cxxopts::value<int>(norder)->default_value("4"))
    ("nodd", "Number of odd basis functions each harmonic subspace",
     cxxopts::value<int>(nodd)->default_value("-1"))
    ("NUMX", "number of radial table entries",
     cxxopts::value<int>(numx)->default_value("128"))
    ("NUMY", "number of vertical table entries",
     cxxopts::value<int>(numy)->default_value("64"))
    ("outr", "number of radial points for output",
     cxxopts::value<int>(OUTR)->default_value("40"))
    ("outz", "number of vertical points for output",
     cxxopts::value<int>(OUTZ)->default_value("40"))
    ("cachefile","cachefile name",
     cxxopts::value<std::string>(CACHEFILE)->default_value(".eof.cache.file"))
    ("outtag", "outtag for basis files",
     cxxopts::value<std::string>(OUTTAG)->default_value("basis"))
    ("density", "compute density",
     cxxopts::value<bool>(DENS)->default_value("true"))
    ("logscale", "logscale for output basis",
     cxxopts::value<bool>(LOGSC)->default_value("false"))
    ;
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0) std::cout << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }

  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAPR       = 1;
  EmpCylSL::CMAPZ       = 1;
  EmpCylSL::logarithmic = true;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::NOUT        = norder;

				// Create expansion
				//
  EmpCylSL ortho(nmax, lmax, mmax, norder, rscale, vscale, nodd, CACHEFILE);
    
  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==1) {
    if (myid==0) {
      std::cout << "EOF parameter match, continuing . . ." << std::flush;
      ortho.dump_basis(OUTTAG, 0, RMAX);
      std::cout << " done" << std::endl;
    }
  } else {
    if (myid==0)
      std::cout << "EOF parameter mismatch . . . quitting" << std::endl;
  }
  
  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

