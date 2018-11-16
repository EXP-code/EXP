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

				// BOOST stuff
#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp> 

namespace po = boost::program_options;
namespace pt = boost::property_tree;


                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <Vector.h>
#include <numerical.h>
#include "Particle.h"
#include <PSP.H>
#include <interp.h>
#include <EmpOrth9thd.h>

#include <localmpi.h>
#include <ProgramParam.H>
#include <foarray.H>

#include <VtkGrid.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

typedef boost::shared_ptr<PSPDump> PSPDumpPtr;

const std::string overview = "Compute and print PCA basis for rendering";

				// Variables not used but needed for linking
int VERBOSE = 4;
int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;
  
				// Globals

int
main(int argc, char **argv)
{
  double RMAX, ZMAX;
  int    OUTR, OUTZ, lmax, mmax, nmax, norder, numx, numy;
  std::string CACHEFILE, OUTTAG;
  double rcylmin, rcylmax, rscale, vscale;
  bool   DENS, LOGSC;

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("RMAX,R",
     po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("ZMAX,Z",
     po::value<double>(&ZMAX)->default_value(0.01),
     "maximum height for output")
    ("rcylmin",
     po::value<double>(&rcylmin)->default_value(0.001),
     "minimum radius for cylindrical basis table")
    ("rcylmax",
     po::value<double>(&rcylmax)->default_value(20.0),
     "maximum radius for cylindrical basis table")
    ("rscale",
     po::value<double>(&rscale)->default_value(0.01), 
     "radial scale length for basis expansion")
    ("vscale",
     po::value<double>(&vscale)->default_value(0.001), 
     "vertical scale length for basis expansion")
    ("lmax",
     po::value<int>(&lmax)->default_value(36), 
     "maximum harmonic order for spherical expansion")
    ("nmax",
     po::value<int>(&nmax)->default_value(8),
     "maximum harmonic order for spherical expansion")
    ("mmax",
     po::value<int>(&mmax)->default_value(4), 
     "maximum azimuthal harmonic order for cylindrical expansion")
    ("norder",
     po::value<int>(&norder)->default_value(4), 
     "maximum radial order for each harmonic subspace")
    ("NUMX",
     po::value<int>(&numx)->default_value(128), 
     "number of radial table entries")
    ("NUMY",
     po::value<int>(&numy)->default_value(64), 
     "number of vertical table entries")
    ("outr",
     po::value<int>(&OUTR)->default_value(40), 
     "number of radial points for output")
    ("outz",
     po::value<int>(&OUTZ)->default_value(40), 
     "number of vertical points for output")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("outtag",
     po::value<std::string>(&OUTTAG)->default_value("basis"),
     "outtag for basis files")
    ("density",
     po::value<bool>(&DENS)->default_value(true),
     "compute density")
    ("logscale",
     po::value<bool>(&LOGSC)->default_value(false),
     "logscale for output basis")
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << overview << std::endl;
    std::cout << desc     << std::endl;
    return 1;
  }

  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAP        = true;
  EmpCylSL::logarithmic = true;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::NOUT        = norder;
  EmpCylSL::CACHEFILE   = CACHEFILE;

				// Create expansion
				//
  EmpCylSL ortho(nmax, lmax, mmax, norder, rscale, vscale);
    
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==1) {
    if (myid==0) {
      std::cout << "EOF parameter match, continuing" << std::flush;
      ortho.dump_basis(OUTTAG, 0, RMAX);
      std::cout << "done" << std::endl;
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

