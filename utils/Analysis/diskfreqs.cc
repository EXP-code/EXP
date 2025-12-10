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
#include <memory>
#include <cmath>
#include <string>

                                // System libs
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include "numerical.H"
#include "interp.H"
#include "EmpCylSL.H"
#include "SphSL.H"

				// Local coefficient classes
#include "Coefs.H"

#include "localmpi.H"
#include "foarray.H"
#include "cxxopts.H"
#include "libvars.H"

const std::string overview = "Compute azimuthal and vertical disk frequencies from coefficients";

int
main(int argc, char **argv)
{
  int numx, numy, lmax=36, mmax, nmax, norder, cmapr, cmapz, numr, nodd=-1;
  double rcylmin, rcylmax, rscale, vscale, RMAX, Tmin, Tmax, dT, H, eps;
  bool verbose = false, mask = false, ignore, logl;
  std::string CACHEFILE, COEFFILE, COEFFILE2, MODEL, OUTFILE, fileType, filePrefix;

  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);

  // ==================================================
  // Parse Command line
  // ==================================================

  cxxopts::Options options("diskfreqs", overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v,verbose", "verbose output")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("P,prefix", "prefix for phase-space files",
     cxxopts::value<std::string>(filePrefix)->default_value("OUT"))
    ("t,Tmin", "Minimum time point for freq evaluation",
     cxxopts::value<double>(Tmin)->default_value("0.0"))
    ("T,Tmax", "Maximum time point for freq evaluation",
     cxxopts::value<double>(Tmax)->default_value("6.0"))
    ("d,dT", "Delta time point for freq evaluation",
     cxxopts::value<double>(dT)->default_value("0.1"))
    ("R,RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("0.1"))
    ("N,numr", "number of output radial grid points",
     cxxopts::value<int>(numr)->default_value("40"))
    ("rscale", "radial scale length for basis expansion",
     cxxopts::value<double>(rscale)->default_value("0.01"))
    ("vscale", "vertical scale length for basis expansion",
     cxxopts::value<double>(vscale)->default_value("0.001"))
    ("H,H", "disk scale height for computing numerical derivative of fz",
     cxxopts::value<double>(H)->default_value("0.002"))
    ("eps", "disk scale height fraction for computing numerical derivative of fz",
     cxxopts::value<double>(eps)->default_value("0.1"))
    ("cachefile", "cachefile name",
     cxxopts::value<std::string>(CACHEFILE)->default_value(".eof.cache.file"))
    ("coeffile", "Disk coefficient file name",
     cxxopts::value<std::string>(COEFFILE))
    ("coeffile2", "Halo coefficient file name",
     cxxopts::value<std::string>(COEFFILE2))
    ("model", "Spherical basis model file",
     cxxopts::value<std::string>(MODEL))
    ("logl", "use logarithmic radius scale in cylindrical grid computation",
     cxxopts::value<bool>(logl)->default_value("true"))
    ("ignore", "rebuild EOF grid if input parameters do not match the cachefile",
     cxxopts::value<bool>(ignore)->default_value("false"))
    ("o,output", "output data table",
     cxxopts::value<std::string>(OUTFILE)->default_value("diskfreqs"))
    ;

  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << std::string(60, '-') << std::endl;
    std::cout << options.help() << std::endl;
    std::cout << std::string(60, '-') << std::endl << std::endl;
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
      auto buf = std::make_unique<char[]>(ssize+1);
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
      if (node["nodd"])
	nodd  = node["nodd"  ].as<int>();
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
      in.read((char *)&tmp,     sizeof(int)); 
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

				// Create expansion
				//
  EmpCylSL ortho_disk(nmax, lmax, mmax, norder, rscale, vscale, nodd, CACHEFILE);
    
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

  auto halo = std::make_shared<SphericalModelTable>(MODEL);
  SphSL::NUMR = 4000;
  int LMAX = 1, NMAX = coefsH.begin()->second->coefs.cols();
  SphSL ortho_halo(halo, LMAX, NMAX);

    
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

    Eigen::MatrixXd sphcoef((LMAX+1)*(LMAX+1), NMAX);
    sphcoef.setZero();		// Need this?

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
      for (int n=0; n<NMAX; n++) sphcoef(0, n) = itH->coefs(0, n);
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

