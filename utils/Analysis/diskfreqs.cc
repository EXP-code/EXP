/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute disk frequencies
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
 *  MDW 09/01/20
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
#include <boost/make_unique.hpp>
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
#include <EmpCylSL.h>
#include <SphereSL.H>

				// Local coefficient classes
#include "Coefs.H"

#include <localmpi.h>
#include <foarray.H>

const std::string overview = "Compute azimuthal and vertical disk frequencies from coefficients";

				// Variables not used but needed for linking
int VERBOSE = 4;
int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
std::vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
std::string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;
  

int
main(int argc, char **argv)
{
  int nice, numx, numy, lmax, mmax, nmax, norder, LMAX, NMAX;
  int initc, partc, beg, end, stride, init, cmapr, cmapz;
  double rcylmin, rcylmax, rscale, vscale, RMAX;
  bool DENS, PCA, PVD, verbose = false, mask = false, ignore, logl;
  std::string CACHEFILE, COEFFILE, COEFFILE2, MODEL;

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("mask,b",
     "blank empty cells")
    ("nice",
     po::value<int>(&nice)->default_value(0), 
     "number of bins in x direction")
    ("RMAX,R",
     po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("rcylmin",
     po::value<double>(&rcylmin)->default_value(0.001),
     "minimum radius for cylindrical basis table")
    ("rcylmax",
     po::value<double>(&rcylmax)->default_value(20.0),
     "maximum radius for cylindrical basis table")
    ("NUMX",
     po::value<int>(&numx)->default_value(128), 
     "number of radial table entries")
    ("NUMY",
     po::value<int>(&numy)->default_value(64), 
     "number of vertical table entries")
    ("rscale",
     po::value<double>(&rscale)->default_value(0.01), 
     "radial scale length for basis expansion")
    ("vscale",
     po::value<double>(&vscale)->default_value(0.001), 
     "vertical scale length for basis expansion")
    ("LMAX",
     po::value<int>(&LMAX)->default_value(1), 
     "maximum harmonic order for spherical basis expansion")
    ("NMAX",
     po::value<int>(&NMAX)->default_value(8),
     "maximum radial order for spherical basis expansion")
    ("lmax",
     po::value<int>(&lmax)->default_value(36), 
     "maximum harmonic order for spherical expansion for disk basis")
    ("nmax",
     po::value<int>(&nmax)->default_value(8),
     "maximum radial order for spherical expansion for disk basis")
    ("mmax",
     po::value<int>(&mmax)->default_value(4), 
     "maximum azimuthal harmonic order for cylindrical expansion")
    ("norder",
     po::value<int>(&norder)->default_value(4), 
     "maximum radial order for each harmonic subspace")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("coeffile",
     po::value<std::string>(&COEFFILE),
     "Disk coefficient file name")
    ("coeffile2",
     po::value<std::string>(&COEFFILE2),
     "Halo coefficient file name")
    ("cmapr",
     po::value<int>(&cmapr)->default_value(1),
     "Radial coordinate mapping type for cylindrical grid (0=none, 1=rational fct)")
    ("cmapz",
     po::value<int>(&cmapz)->default_value(1),
     "Vertical coordinate mapping type for cylindrical grid (0=none, 1=sech, 2=power in z")
    ("model",
     po::value<std::string>(&MODEL),
     "Spherical basis model file")
    ("logl",
     po::value<bool>(&logl)->default_value(true),
     "use logarithmic radius scale in cylindrical grid computation")
    ("ignore",
     po::value<bool>(&ignore)->default_value(false),
     "rebuild EOF grid if input parameters do not match the cachefile")
    ("runtag",
     po::value<std::string>(&runtag)->default_value("run1"),
     "runtag for phase space files")
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
    std::cout << std::string(60, '-') << std::endl;
    std::cout << overview << std::endl;
    std::cout << std::string(60, '-') << std::endl << std::endl;
    std::cout << desc     << std::endl;
    return 1;
  }
 
  if (vm.count("verbose")) verbose = true;

  if (vm.count("mask")) mask = true;


#ifdef DEBUG
  sleep(20);
#endif  

  // Okay here is the plan
  // ---------------------
  // (1) Read in and initialize the disk coefficients and basis
  // (2) Read in and initialize the halo coefficients and basis
  // (3) For output, we basically want a matrix of Omega(R, t)
  // (4) For each desired time (perhaps have the user pick a dT),
  //     then simply do the evalutions of Omega_R = \sqrt{-f_R/R}
  //     and Omega_z = \sqrt{d^2V/dz^2}.

  std::ifstream in(CACHEFILE);
  if (!in) {
    std::cerr << "Error opening cachefile named <" 
	      << CACHEFILE << ">, quitting" << std::endl;
  } else {

    // Attempt to read magic number
    //
    unsigned int tmagic;
    in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

    //! Basis magic number
    const unsigned int hmagic = 0xc0a57a1;
      
    if (tmagic == hmagic) {
      // YAML size
      //
      unsigned ssize;
      in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));
	
      // Make and read char buffer
      //
      auto buf = boost::make_unique<char[]>(ssize+1);
      in.read(buf.get(), ssize);
      buf[ssize] = 0;		// Null terminate
      
      YAML::Node node;
      
      try {
	node = YAML::Load(buf.get());
      }
      catch (YAML::Exception& error) {
	if (myid)
	  std::cerr << "YAML: error parsing <" << buf.get() << "> "
		    << "in " << __FILE__ << ":" << __LINE__ << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
      }

      // Get parameters
      //
      mmax    = node["mmax"  ].as<int>();
      numx    = node["numx"  ].as<int>();
      numy    = node["numy"  ].as<int>();
      nmax    = node["nmax"  ].as<int>();
      norder  = node["norder"].as<int>();
      DENS    = node["dens"  ].as<bool>();
      if (node["cmap"])
	cmapr = node["cmap"  ].as<int>();
      else
	cmapr = node["cmapr" ].as<int>();
      if (node["cmapz"])
	cmapz = node["cmapz"  ].as<int>();
      rcylmin = node["rmin"  ].as<double>();
      rcylmax = node["rmax"  ].as<double>();
      rscale  = node["ascl"  ].as<double>();
      vscale  = node["hscl"  ].as<double>();
      
    } else {
				// Rewind file
      in.clear();
      in.seekg(0);
      
      int tmp;
      
      in.read((char *)&mmax,    sizeof(int));
      in.read((char *)&numx,    sizeof(int));
      in.read((char *)&numy,    sizeof(int));
      in.read((char *)&nmax,    sizeof(int));
      in.read((char *)&norder,  sizeof(int));
      in.read((char *)&DENS,    sizeof(int)); 
      in.read((char *)&cmapr,   sizeof(int)); 
      in.read((char *)&rcylmin, sizeof(double));
      in.read((char *)&rcylmax, sizeof(double));
      in.read((char *)&rscale,  sizeof(double));
      in.read((char *)&vscale,  sizeof(double));
    }
  }
  
  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAPR       = cmapr;
  EmpCylSL::CMAPZ       = cmapz;
  EmpCylSL::logarithmic = logl;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::CACHEFILE   = CACHEFILE;

				// Create expansion
				//
  EmpCylSL ortho_disk(nmax, lmax, mmax, norder, rscale, vscale);
    
  // ==================================================
  // Initialize disk create basis
  // ==================================================
  
  if (ortho_disk.read_cache()==0) {
    std::cerr << "Error eading cache file <" << CACHEFILE << ">"
	      << std::endl;
    exit(-3);
  }
  

  // ==================================================
  // Open disk input coefficient file
  // ==================================================
  
  std::ifstream coefs_disk(COEFFILE);
  if (not coefs_disk) {
    std::cerr << "Could not open coefficient file <" << COEFFILE
	      << "> . . . quitting" << std::endl;
    exit(-4);
  }

  std::map<double, CoefPtr> coefsD;

  while (coefs_disk) {
    CoefPtr c = std::make_shared<Coefs>();
    if (not c->read(coefs_disk, verbose)) break;

    coefsD[c->time] = c;
  }
  
  // ==================================================
  // Initialize halo basis
  // ==================================================
  
  // ==================================================
  // Make SL expansion
  // ==================================================

  SphericalModelTable halo(MODEL);
  SphereSL::NUMR = 4000;
  SphereSL ortho_halo(&halo, LMAX, NMAX);


  // ==================================================
  // Open disk input coefficient file
  // ==================================================
  
  std::ifstream coefs_halo(COEFFILE2);
  if (not coefs_halo) {
    std::cerr << "Could not open coefficient file <" << COEFFILE
	      << "> . . . quitting" << std::endl;
    exit(-5);
  }

  std::map<double, SphCoefsPtr> coefsH;

  while (in) {
    SphCoefsPtr c = std::make_shared<SphCoefs>();
    if (not c->read(coefs_halo, verbose)) break;

    coefsH[c->header.tnow] = c;
  }
    
  // DONE
  //
  return 0;
}

