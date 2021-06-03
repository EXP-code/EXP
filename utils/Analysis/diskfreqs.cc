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
#include <algorithm>
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
  int numx, numy, lmax=36, mmax, nmax, norder, cmapr, cmapz, numr;
  double rcylmin, rcylmax, rscale, vscale, RMAX, Tmin, Tmax, dT, H, eps;
  bool DENS, verbose = false, mask = false, ignore, logl;
  std::string CACHEFILE, COEFFILE, COEFFILE2, MODEL, OUTFILE;

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("Tmin,t",
     po::value<double>(&Tmin)->default_value(0.0),
     "Minimum time point for freq evaluation")
    ("Tmax,T",
     po::value<double>(&Tmax)->default_value(6.0),
     "Maximum time point for freq evaluation")
    ("dT,d",
     po::value<double>(&dT)->default_value(0.1),
     "Delta time point for freq evaluation")
    ("RMAX,R",
     po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("numr,N",
     po::value<int>(&numr)->default_value(40),
     "number of output radial grid points")
    ("rscale",
     po::value<double>(&rscale)->default_value(0.01), 
     "radial scale length for basis expansion")
    ("vscale",
     po::value<double>(&vscale)->default_value(0.001), 
     "vertical scale length for basis expansion")
    ("H,H",
     po::value<double>(&H)->default_value(0.002), 
     "disk scale height for computing numerical derivative of fz")
    ("eps",
     po::value<double>(&eps)->default_value(0.1), 
     "disk scale height fraction for computing numerical derivative of fz")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("coeffile",
     po::value<std::string>(&COEFFILE),
     "Disk coefficient file name")
    ("coeffile2",
     po::value<std::string>(&COEFFILE2),
     "Halo coefficient file name")
    ("model",
     po::value<std::string>(&MODEL),
     "Spherical basis model file")
    ("logl",
     po::value<bool>(&logl)->default_value(true),
     "use logarithmic radius scale in cylindrical grid computation")
    ("ignore",
     po::value<bool>(&ignore)->default_value(false),
     "rebuild EOF grid if input parameters do not match the cachefile")
    ("output,o",
     po::value<std::string>(&OUTFILE)->default_value("diskfreqs"),
     "output data table")
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

  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
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
    MPI_Finalize();
    exit(-3);
  }
  

  // ==================================================
  // Open disk input coefficient file
  // ==================================================
  
  std::ifstream coefs_disk(COEFFILE);
  if (not coefs_disk) {
    std::cerr << "Could not open coefficient file <" << COEFFILE
	      << "> . . . quitting" << std::endl;
    MPI_Finalize();
    exit(-4);
  }

  std::map<double, CylCoefsPtr> coefsD;

  while (coefs_disk) {
    CylCoefsPtr c = std::make_shared<CylCoefs>();
    if (not c->read(coefs_disk, verbose)) break;

    coefsD[c->time] = c;
  }
  
  // ==================================================
  // Open disk input coefficient file
  // ==================================================
  
  std::ifstream coefs_halo(COEFFILE2);
  if (not coefs_halo) {
    std::cerr << "Could not open coefficient file <" << COEFFILE
	      << "> . . . quitting" << std::endl;
    MPI_Finalize();
    exit(-5);
  }

  std::map<double, SphCoefsPtr> coefsH;

  while (in) {
    SphCoefsPtr c = std::make_shared<SphCoefs>();
    if (not c->read(coefs_halo, verbose)) break;

    coefsH[c->header.tnow] = c;
  }

  // ==================================================
  // Make SL expansion
  // ==================================================

  SphericalModelTable halo(MODEL);
  SphereSL::NUMR = 4000;
  int LMAX = 1, NMAX = coefsH.begin()->second->coefs[0].size();
  SphereSL ortho_halo(&halo, LMAX, NMAX);

    
  // Begin grid creation
  //
  std::vector<double> times, rgrid(numr);

  for (int i=0; i<numr; i++) rgrid[i] = RMAX*(0.5+i)/numr;

  // Create time grid
  //
  double tcur = Tmin;
  while (tcur < Tmax) {
    auto itD = coefsD.lower_bound(tcur);
    auto itH = coefsH.lower_bound(tcur);
    if (fabs(itD->first - itH->first) > 1.0e-8) {
      std::cerr << "Time mismatch for T=" << tcur
		<< ", Disk=" << itD->first << ", Halo=" << itH->first
		<< std::endl;
    } else {
      if (times.size()==0 or (itD->first - times.back() > 1.0e-8))
	times.push_back(itD->first);
    }
    tcur += dT;
  }

  // Open output file
  //
  std::ofstream outR(OUTFILE + ".omp");
  std::ofstream outZ(OUTFILE + ".omz");

  if (outR and outZ) {

    Matrix sphcoef(0, LMAX*(LMAX+2), 1, NMAX);
    sphcoef.zero();		// Need this?

    // Radial grid points
    //
    for (auto r : rgrid) {
      outR << std::setw(16) << r;
      outZ << std::setw(16) << r;
    }
    outR << std::endl;
    outZ << std::endl;

    // Compute values
    //
    for (auto t : times) {

      outR << std::setw(16) << t;
      outZ << std::setw(16) << t;

      // Set disk coefficients for m=0 only
      //
      auto itD = coefsD.find(t)->second;
      
      ortho_disk.set_coefs(0, itD->cos_c[0], itD->sin_c[0], true);

      // Set halo coefficients for l=m=0 only
      //
      auto itH = coefsH.find(t)->second;
      for (int n=0; n<NMAX; n++) sphcoef[0][n+1] = itH->coefs[0][n];
      ortho_halo.install_coefs(sphcoef);

      for (auto r : rgrid) {
	// Now compute forces
	double z = 0.0, phi = 0.0, rforce = 0.0, costh = 0.0;

	double p0, p1, fr, fz, fp, d0, d1, ft;
	ortho_disk.accumulated_eval(r, z, phi, p0, p1, fr, fz, fp);

	rforce += fr;

	ortho_halo.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);

	rforce -= fr;

	if (rforce<0.0) rforce = sqrt(-rforce/r);
	else rforce = 0.0;

	outR << std::setw(16) << rforce;

	double dfz = 0.0;
	ortho_disk.accumulated_eval(r, z+eps*H, phi, p0, p1, fr, fz, fp);
	dfz += fz;
	ortho_disk.accumulated_eval(r, z-eps*H, phi, p0, p1, fr, fz, fp);
	dfz -= fz;
	dfz /= 2.0*eps*H;
	
	if (dfz<0.0) dfz = sqrt(-dfz);
	else dfz = 0.0;

	outZ << std::setw(16) << dfz;

      }
      outR << std::endl;
      outZ << std::endl;
    }
  }

  MPI_Finalize();

  // DONE
  //
  return 0;
}

