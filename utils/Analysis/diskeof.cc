/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Create EOF from a sequence of PSP files
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
 *  MDW 12/11/20
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

#include <config.h>

#ifdef HAVE_LIBPNGPP
#include <ColorGradient.H>	// For PNG images
#endif

#include <yaml-cpp/yaml.h>	// YAML support

#include <Eigen/Eigen>		// Eigen 3

namespace po = boost::program_options;
namespace pt = boost::property_tree;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <Vector.h>
#include <numerical.h>
#include "Particle.h"
#include <PSP2.H>
#include <interp.h>
#include <EmpCylSL.h>

#include <localmpi.h>
#include <foarray.H>

#include <VtkGrid.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

using CoefArray  = std::vector<Eigen::VectorXd>;
using DvarArray  = std::vector<Eigen::MatrixXd>;

const std::string overview = "Compute new EOF basis PSP phase-space output files";

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
  int lmax=64, mmax, nmax, norder, numx, numy, cmapr=1, cmapz=1;
  double rcylmin, rcylmax, rscale, vscale;
  std::string CACHEFILE, COEFFILE, cname;
  int beg, end, stride;

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("PNG",
     "PNG matrix output")
    ("beg",
     po::value<int>(&beg)->default_value(0),
     "initial PSP index")
    ("end",
     po::value<int>(&end)->default_value(99999),
     "final PSP index")
    ("stride",
     po::value<int>(&stride)->default_value(1),
     "PSP index stride")
    ("outdir",
     po::value<std::string>(&outdir)->default_value("."),
     "Output directory path")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("coeffile",
     po::value<std::string>(&COEFFILE),
     "coefficient output file name")
    ("cname",
     po::value<std::string>(&cname)->default_value("star"),
     "component name")
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
 
  bool PNG = false;
  if (vm.count("PNG")) PNG = true;

#ifdef DEBUG
  sleep(20);
#endif  

  bool DENS = false;
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  // ==================================================
  // Read basis cache
  // ==================================================

  std::ifstream in(CACHEFILE);
  if (!in) {
    if (myid==0)
      std::cerr << "Error opening cachefile named <" 
		<< CACHEFILE << "> . . . I quit!"
		<< std::endl;
    MPI_Finalize();
    exit(-1);

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
      in.read((char *)&tmp,     sizeof(int));   
      in.read((char *)&cmapr,   sizeof(int)); 
      in.read((char *)&rcylmin, sizeof(double));
      in.read((char *)&rcylmax, sizeof(double));
      in.read((char *)&rscale,  sizeof(double));
      in.read((char *)&vscale,  sizeof(double));

      if (tmp) DENS = true;
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

  // Create expansion instance
  //
  EmpCylSL ortho(nmax, lmax, mmax, norder, rscale, vscale);
    
  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==0) {
    if (myid==0)
      std::cerr << "Error creating disk basis using <"
		<< CACHEFILE << "> . . . I quit!"
		<< std::endl;
    MPI_Finalize();
    exit(-1);
  }
  
  // ==================================================
  // Create new basis and coefficients
  // ==================================================
  
  // First create storage
  //
  DvarArray D(mmax+1);
  for (auto & mat : D) {
    mat = Eigen::MatrixXd::Zero(norder, norder);
  }

  std::vector<CoefArray> coefsC, coefsS;
  std::vector<double>    times;
  std::vector<int>       indices;

  for (int indx=beg; indx<=end; indx+=stride) {

    // ==================================================
    // PSP input stream
    // ==================================================

    int iok = 1;
    std::ostringstream s1;
    s1 << "OUT." << runtag << "."
       << std::setw(5) << std::setfill('0') << indx;
      
    if (myid==0) {
      std::ifstream in1(s1.str());	// Now, try to open a new one . . . 
      if (!in1) {
	cerr << "Error opening <" << s1.str() << ">" << endl;
	iok = 0;
      }
    }
    
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (iok==0) break;
    
    // ==================================================
    // Open frame list
    // ==================================================
    
    auto psp = boost::make_shared<PSPout>(s1.str());
    
    if (psp->GetNamed(cname) == 0) {
      if (myid==0) std::cout << "Could not find component <" << cname
			     << "> in PSP file <" << s1.str() << ">"
			     << std::endl;
      MPI_Finalize();
      exit(-2);
    }

    tnow = psp->CurrentTime();

    if (myid==0) {
      cout << "Beginning disk partition [time=" << tnow
	   << ", index=" << indx << "] . . . "  << flush;
    }

    times.push_back(tnow);
    indices.push_back(indx);
    
    SParticle *p = psp->GetParticle();

    CoefArray coefC(mmax + 1);	// Per snapshot storage for
    CoefArray coefS(mmax + 1);	// coefficients

    for (auto & v : coefC) v = Eigen::VectorXd::Zero(norder);
    for (auto & v : coefS) v = Eigen::VectorXd::Zero(norder);
    
    unsigned count = 0;
    
    EmpCylSL::ContribArray retC, retS;

    while (p) {
      if (count++ % numprocs == myid) {
	// Only need mass and position
	//
	double m = p->mass();
	double x = p->pos(0);
	double y = p->pos(1);
	double z = p->pos(2);
	double p = atan2(y, x);

	// Get coefficient contribution for this particle
	//
	ortho.getPotParticle(x, y, z, retC, retS);
	
	// Accumulate coefficients and D matrix
	//
	for (int mm=0; mm<=mmax; mm++) {
	  for (int n1=0; n1<norder; n1++) {
	    // Coefficient contribution
	    coefC[mm](n1) += m * retC[mm](n1);
	    if (mm) coefS[mm](n1) += m * retS[mm](n1);

	    // Modulus for index n1
	    double mod1 = retC[mm](n1) * retC[mm](n1);
	    if (mm) mod1 += retS[mm](n1) * retS[mm](n1);

	    for (int n2=0; n2<norder; n2++) {
	      // Modulus for index n2
	      double mod2 = retC[mm](n2) * retC[mm](n2);
	      if (mm) mod2 += retS[mm](n2) * retS[mm](n2);

	      D[mm](n1, n2) += m * sqrt(mod1 * mod2);
	    }
	  }
	}
      }
      p = psp->NextParticle();
    }

    coefsC.push_back(coefC);
    coefsS.push_back(coefS);

    if (myid==0) std::cout << "done" << std::endl;

  } // PSP loop

  // Full reduction
  //
  for (auto & mat : D) {
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), mat.size(), MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
  }

  for (auto & coef : coefsC) {
    for (auto & d : coef) 
      MPI_Allreduce(MPI_IN_PLACE, d.data(), d.size(), MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);
  }

  for (auto & coef : coefsS) {
    for (auto & d : coef) 
      MPI_Allreduce(MPI_IN_PLACE, d.data(), d.size(), MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);
  }

  // Okay, now compute SVD
  //
  std::vector<Eigen::VectorXd> S;
  std::vector<Eigen::MatrixXd> U;

  for (int mm=0; mm<=mmax; mm++) {

    Eigen::JacobiSVD<Eigen::MatrixXd>
      svd(D[mm], Eigen::ComputeThinU | Eigen::ComputeThinV);

    if (myid==0)
      std::cout << "Singular values for m=" << mm << std::endl
		<< svd.singularValues() << std::endl;

    S.push_back(svd.singularValues());
    auto u = svd.matrixU();
    U.push_back(u);


#ifdef HAVE_LIBPNGPP
    if (PNG and myid==0) {

      const int minSize = 600;
      int ndup = 1;

      if (norder < minSize) ndup = minSize/norder   + 1;

      png::image< png::rgb_pixel > image(norder*ndup, norder*ndup);
      ColorGradient color;
      // color.createFiveColorHeatMapGradient();
      color.createGrayGradient();
      
      std::ostringstream sout;
      sout << "diskeof_disp." << mm;

      double minV = std::numeric_limits<double>::max();
      double maxV = std::numeric_limits<double>::min();

      for (int n1=0; n1<norder; n1++) {
	for (int n2=0; n2<norder; n2++) {
	  minV = std::min<double>(minV, D[mm](n1, n2));
	  maxV = std::max<double>(maxV, D[mm](n1, n2));
	}
      }

      for (int n1=0; n1<norder; n1++) {
	for (int n2=0; n2<norder; n2++) {
	  png::rgb_pixel cval = color( (D[mm](n1, n2) - minV)/(maxV - minV) );
	  for (size_t xx = n1*ndup; xx < (n1+1)*ndup; xx++) {
	    for (size_t yy = n2*ndup; yy < (n2+1)*ndup; yy++) {
	      image[xx][yy] = cval;
	    }
	  }
	}
      }
      
      image.write(sout.str() + ".png");
      
      sout.str("");
      sout << "diskeof_ev." << mm;

      minV = -1.0;
      maxV =  1.0;

      for (int n1=0; n1<norder; n1++) {
	for (int n2=0; n2<norder; n2++) {
	  png::rgb_pixel cval = color( (u(n1, n2) - minV)/(maxV - minV) );
	    for (size_t xx = n1*ndup; xx < (n1+1)*ndup; xx++) {
	      for (size_t yy = n2*ndup; yy < (n2+1)*ndup; yy++) {
		image[xx][yy] = cval;
	    }
	  }
	}
      }
      
      image.write(sout.str() + ".png");
    }
#endif
  }

  // Make a coefficient file for rotating coefficients in the same
  // format as readcoefs
  //
  if (myid==0) {
    std::ofstream out("diskeofs.coefs");
    std::ofstream org("diskeofs.coefs_orig");
    int ntimes = times.size();
    for (int t=0; t<ntimes; t++) {
      for (int mm=0; mm<=mmax; mm++) {
	out << std::setw(18) << times[t] << std::setw(5) << mm;
	org << std::setw(18) << times[t] << std::setw(5) << mm;

	auto newC = U[mm].transpose() * coefsC[t][mm];
	Eigen::VectorXd newS = Eigen::VectorXd::Zero(norder);
	if (mm) newS = U[mm].transpose() * coefsS[t][mm];

	for (int nn=0; nn<norder; nn++) {
	  out << std::setw(18) << sqrt(newC[nn]*newC[nn] + newS[nn]*newS[nn]) / ntimes;
	}
	out << std::endl;

	for (int nn=0; nn<norder; nn++) {
	  double res = coefsC[t][mm][nn] * coefsC[t][mm][nn] / ntimes;
	  if (mm) res += coefsS[t][mm][nn] * coefsS[t][mm][nn] / ntimes;
	  org << std::setw(18) << sqrt(res);
	}
	org << std::endl;
      }
    }
  }

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

