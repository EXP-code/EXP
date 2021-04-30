/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Kullback-Leibler analysis for cylinder
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
 *  MDW 11/27/20
 *
 ***************************************************************************/

				// C++/STL headers
#include <numeric>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <queue>
#include <map>

				// Eigen3
#include <Eigen/Eigen>

				// Boost stuff

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/make_unique.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <Progress.H>

namespace po = boost::program_options;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <Vector.h>
#include <numerical.h>
#include "Particle.h"
#include <PSP2.H>
#include <interp.h>
#include <massmodel.h>
#include <EmpCylSL.h>
#include <foarray.H>
#include <KDtree.H>

#include <localmpi.h>

#include <yaml-cpp/yaml.h>	// YAML support

// Variables not used but needed for linking
//
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
std::string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;
  
// Helper class
//

class CoefStruct
{

private:
  int mmax, nmax;

public:
  std::vector<std::vector<double>> coefC, coefS;

  CoefStruct(int mmax, int nmax) : mmax(mmax), nmax(nmax)
  {
    coefC.resize(mmax+1);
    coefS.resize(mmax+1);
    for (int m=0; m<=mmax; m++) {
      coefC[m].resize(nmax);
      if (m) coefS[m].resize(nmax);
    }
  }

  void sync(double norm)
  {
    // Sum the norm from all processes
    //
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1,
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Sanity check
    //
    if (norm <= 0.0) norm = 1.0;

    // Apply the normalization
    //
    for (int m=0; m<=mmax; m++) {

      for (auto & v : coefC[m]) v /= norm;

      if (m) {
	for (auto & v : coefS[m]) v /= norm;
      }
    }
  }

};
  

int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, rscale, minSNR0, Hexp;
  int NICE, LMAX, NMAX, NSNR, indx, nbunch, Ndens;
  std::string CACHEFILE, dir("./"), cname, prefix;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Kullback-Leibler analysis for cylindrical models" << std::endl
       << std::string(60, '-') << std::endl << std::endl
       << "Allowed options";
  
  po::options_description desc(sout.str());
  desc.add_options()
    ("help,h",
     "Print this help message")
    ("verbose,v",
     "Verbose and diagnostic output for covariance computation")
    ("truncate,t",
     "Use Truncate method for SNR trimming rather than the default Hall")
    ("debug,",
     "Debug max values")
    ("OUT",
     "assume original, single binary PSP files as input")
    ("SPL",
     "assume new split binary PSP files as input")
    ("LOG",
     "log scaling for SNR")
    ("Hall",
     "use Hall smoothing for SNR trim")
    ("Ndens,K",             po::value<int>(&Ndens)->default_value(32),
     "KD density estimate count (use 0 for expansion estimate)")
    ("NICE",                po::value<int>(&NICE)->default_value(0),
     "system priority")
    ("LMAX",                po::value<int>(&LMAX)->default_value(36),
     "Maximum harmonic order for spherical expansion")
    ("NSNR, N",             po::value<int>(&NSNR)->default_value(20),
     "Number of SNR evaluations")
    ("minSNR",              po::value<double>(&minSNR0),
     "minimum SNR value for loop output")
    ("Hexp",                po::value<double>(&Hexp)->default_value(1.0),
     "default Hall smoothing exponent")
    ("prefix",              po::value<string>(&prefix)->default_value("crossval"),
     "Filename prefix")
    ("runtag",              po::value<string>(&runtag)->default_value("run1"),
     "Phase space file")
    ("outdir",              po::value<string>(&outdir)->default_value("."),
     "Output directory path")
    ("indx",                po::value<int>(&indx)->default_value(0),
     "PSP index")
    ("nbunch",              po::value<int>(&nbunch)->default_value(-1),
     "Desired bunch size (default: sqrt(nbod) if value is < 0)")
    ("dir,d",               po::value<std::string>(&dir),
     "directory for SPL files")
    ("ignore",
     po::value<bool>(&ignore)->default_value(false),
     "rebuild EOF grid if input parameters do not match the cachefile")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("cname",
     po::value<std::string>(&cname)->default_value("star disk"),
     "component name")
    ;
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    if (myid==0) std::cout << std::endl << desc << std::endl;
    MPI_Finalize();
    return 0;
  }

  bool SPL = false;
  if (vm.count("SPL")) SPL = true;
  if (vm.count("OUT")) SPL = false;

  bool LOG = false;
  if (vm.count("LOG")) {
    LOG = true;
  }

  bool Hall = false;
  if (vm.count("Hall")) Hall = true;

  bool verbose = false;
  if (vm.count("verbose")) verbose = true;

  bool debug = false;
  if (vm.count("debug")) debug = true;

  // ==================================================
  // Nice process
  // ==================================================

  if (NICE>0)
    setpriority(PRIO_PROCESS, 0, NICE);

  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  int mmax, numx, numy, norder, cmapr, cmapz;
  double rcylmin, rcylmax, vscale;
  bool DENS;

  if (not ignore) {

    std::ifstream in(CACHEFILE);
    if (!in) {
      std::cerr << "Error opening cachefile named <" 
		<< CACHEFILE << "> . . ."
		<< std::endl
		<< "I will build <" << CACHEFILE
		<< "> but it will take some time."
		<< std::endl
		<< "If this is NOT what you want, "
		<< "stop this routine and specify the correct file."
		<< std::endl;
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
	NMAX    = node["nmax"  ].as<int>();
	norder  = node["norder"].as<int>();
	DENS    = node["dens"  ].as<bool>();
	if (node["cmap"])
	  cmapr = node["cmap"  ].as<int>();
	else
	  cmapr = node["cmapr" ].as<int>();
	if (node["cmapz"])
	  cmapz = node["cmapz" ].as<int>();
	rcylmin = node["rmin"  ].as<double>();
	rcylmax = node["rmax"  ].as<double>();
	rscale  = node["ascl"  ].as<double>();
	vscale  = node["hscl"  ].as<double>();
	
      } else {
				// Rewind file
	in.clear();
	in.seekg(0);

	int idens;
    
	in.read((char *)&mmax,    sizeof(int));
	in.read((char *)&numx,    sizeof(int));
	in.read((char *)&numy,    sizeof(int));
	in.read((char *)&NMAX,    sizeof(int));
	in.read((char *)&norder,  sizeof(int));
	in.read((char *)&idens,   sizeof(int)); 
	in.read((char *)&cmapr,   sizeof(int)); 
	in.read((char *)&rcylmin, sizeof(double));
	in.read((char *)&rcylmax, sizeof(double));
	in.read((char *)&rscale,  sizeof(double));
	in.read((char *)&vscale,  sizeof(double));

	if (idens) DENS = true;
      }
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
  EmpCylSL::PCAVAR      = true;
  EmpCylSL::PCADRY      = true;

				// Create expansion
				//
  EmpCylSL ortho0(NMAX, LMAX, mmax, norder, rscale, vscale);
  EmpCylSL ortho1(NMAX, LMAX, mmax, norder, rscale, vscale);
    
				// Set smoothing type to Truncate or
				// Hall (default)
  EmpCylSL::HEXP = Hexp;
  if (vm.count("truncate")) {
    ortho0.setTK("Truncate");
    ortho1.setTK("Truncate");
  } else {
    ortho0.setTK("Hall");
    ortho1.setTK("Hall");
  }

  PSPptr psp;
  
  if (ortho0.read_cache()==0 or ortho1.read_cache()==0) {
    std::cout << "Could not read cache file <" << CACHEFILE << ">"
	      << " . . . quitting" << std::endl;
    MPI_Finalize();
    exit(0);
  }
  
  // ==================================================
  // Phase space
  // ==================================================

  std::string file;

#ifdef DEBUG
  std::cout << "[" << myid << "] Begin phase -space loop" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif	      

  // ==================================================
  // PSP input stream
  // ==================================================

  int iok = 1;
  std::ostringstream s1;
  
  s1.str("");		// Clear stringstream
  if (SPL) s1 << "SPL.";
  else     s1 << "OUT.";
  s1 << runtag << "."<< std::setw(5) << std::setfill('0') << indx;
      
				// Check for existence of next file
  file = dir + s1.str();
  std::ifstream in(file);
  if (!in) {
    std::cerr << "Error opening <" << file << ">" << endl;
    iok = 0;
  }
  
  if (iok==0) {
    MPI_Finalize();
    exit(-1);
  }

  // ==================================================
  // Open output file
  // ==================================================

  std::ofstream out;
  bool ok = true;
  if (myid==0) {
    out.open(prefix + ".out");
    if (!out) {
      std::cerr << "Error opening output file <" << prefix + ".out" << ">" << std::endl;
      ok = false;
    }
  }
  
  {
    int okay = ok ? 1 : 0;
    MPI_Bcast(&okay, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (okay==0) {
      MPI_Finalize();
      exit(-2);
    }
  }

  // ==================================================
  // Open PSP file
  // ==================================================

  if (SPL) psp = std::make_shared<PSPspl>(s1.str(), dir, true);
  else     psp = std::make_shared<PSPout>(file, true);
  
  tnow = psp->CurrentTime();
  if (myid==0) std::cout << "Beginning partition [time=" << tnow
			 << ", index=" << indx << "] . . . "  << flush;
  
  if (not psp->GetNamed(cname)) {
    if (myid==0) {
      std::cout << "Error finding component named <" << cname << ">" << std::endl;
      psp->PrintSummary(std::cout);
    }
    MPI_Finalize();
    exit(-1);
  }
      
  int nbod = psp->GetNamed(cname)->comp.nbod;

  std::vector<double> KDdens;

  boost::shared_ptr<boost::progress_display> progress;
  if (myid==0) {
    std::cout << std::endl
	      << "Accumulating particle positions . . . "
	      << std::endl;
    progress = boost::make_shared<boost::progress_display>(nbod);
  }

  ortho0.setup_accumulation();
  ortho0.setHall("test", nbod);

  SParticle *p = psp->GetParticle();
  int icnt = 0;
  do {
    if (myid==0) ++(*progress);

    if (icnt++ % numprocs == myid) {
      double R   = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
      double phi = atan2(p->pos(1), p->pos(0));
      ortho0.accumulate(R, p->pos(2), phi, p->mass(), p->indx(), 0, 0, true);
      //                                                         ^  ^  ^
      //                                                         |  |  |
      // Thread id ----------------------------------------------+  |  |
      // Level -----------------------------------------------------+  |
      // Compute covariance -------------------------------------------+
    }
    p = psp->NextParticle();
  } while (p);
  
    
  // This is the kd- NN density estimate; skipped by default for Ndens=0
  //
  if (Ndens) {
    if (myid==0) std::cout << "Computing KD density estimate for " << nbod
			   << " points" << std::endl;

    typedef point <double, 3> point3;
    typedef kdtree<double, 3> tree3;

    std::vector<point3> points;

    // Every node needs to make the tree (a parallel share could be
    // implemented in KDtree.H)
    //
    double KDmass = 0.0;
    for (auto part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
      double ms = part->mass();
      KDmass += ms;
      points.push_back(point3({part->pos(0), part->pos(1), part->pos(2)}, ms));
    }
    
    tree3 tree(points.begin(), points.end());
    
    KDdens.resize(nbod, 0.0);

    int badVol = 0;

    // Share the density computation among the nodes
    //
    for (int k=0; k<points.size(); k++) {
      if (k % numprocs == myid) {
	auto ret = tree.nearestN(points[k], Ndens);
	double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	if (volume>0.0 and KDmass>0.0)
	  KDdens[k] = std::get<1>(ret)/volume/KDmass;
	else badVol++;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, KDdens.data(), nbod,
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << "Finished KD density estimate with " << badVol
		<< " undetermined densities, mass=" << KDmass << std::endl;
      std::cout << "A few densities are: " << KDdens.front()
		<< ", " << KDdens[nbod/2] << ", " << KDdens.back() << std::endl;
    }
  }
  
  //------------------------------------------------------------ 
      
  if (myid==0) cout << "Making coefficients for total . . . " << flush;
  ortho0.make_coefficients(true);
  ortho0.pca_hall(true);
  if (myid==0) std::cout << "done" << endl;
  
  if (myid==0) std::cout << std::endl
			 << "Accumulating particle positions for subsamples . . . "
			 << std::endl;

				// Size of bunch
  int nbunch1 = std::floor(sqrt(nbod));
  if (nbunch>0) nbunch1 = nbod/nbunch; 
				// Number of bunches
  int nbunch0 = nbod/nbunch1;

  double ampfac = 1.0;
  if (nbunch0 > 1) ampfac = 1.0/(nbunch0 - 1);

  p = psp->GetParticle();
  icnt = 0;
    
  std::vector<boost::shared_ptr<CoefStruct>> coefs;

  if (myid==0) {
    progress = boost::make_shared<boost::progress_display>(nbunch0*nbunch1);
  }

  double curMass = 0.0;

  do {
    if (myid==0) ++(*progress);
				// Done processing
    if (icnt >= nbunch0*nbunch1) {
      if (myid==0)
	std::cout << "Finished processing subsamples with n=" << icnt
		  << " N=" << nbod << std::endl;
      break;
    }
      
				// Start a new bunch?
    if (icnt % nbunch1 == 0) {
      if (coefs.size()) {
	ortho1.make_coefficients();
	for (int mm=0; mm<=mmax; mm++) {
	  ortho1.get_coefs(mm,
			   coefs.back()->coefC[mm],
			   coefs.back()->coefS[mm]);
	}
	coefs.back()->sync(curMass);
      }
      coefs.push_back(boost::make_shared<CoefStruct>(mmax, norder));
      ortho1.setup_accumulation();
      curMass = 0.0;
    }
    
				// Particle accumulation is spread
				// between nodes
    if (icnt++ % numprocs == myid) {
      
      double R   = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
      double phi = atan2(p->pos(1), p->pos(0));
      ortho1.accumulate(R, p->pos(2), phi, p->mass(), p->indx(), 0, 0, false);
      //                                                         ^  ^  ^
      //                                                         |  |  |
      // Thread id ----------------------------------------------+  |  |
      // Level -----------------------------------------------------+  |
      // Compute covariance -------------------------------------------+

      curMass += p->mass();	// Accumulate mass per bunch
    }
    p = psp->NextParticle();
  } while (p);
  
  //------------------------------------------------------------ 
      
  if (myid==0) cout << std::endl
		    << "Making coefficients for bunch ["
		    << coefs.size() << "/" << nbunch0 << "]" << flush;

  ortho1.make_coefficients();
  for (int mm=0; mm<=mmax; mm++) {
    ortho1.get_coefs(mm,
		     coefs.back()->coefC[mm],
		     coefs.back()->coefS[mm]);
  }
  coefs.back()->sync(curMass);
  
  if (myid==0) std::cout << "done" << endl;
  
  // This is a debug test
  if (myid==0) {
    std::ofstream test(prefix + ".coeftest");
    for (int n=0; n<norder; n++) {
      test << std::setw(4) << n;
      for (auto & c : coefs)
	test << std::setw(18) << c->coefC[0][n];
      test << std::endl;
    }
  }

  //------------------------------------------------------------ 
    
  if (myid==0) cout << "Beginning SNR loop . . ." << std::endl;


  double minSNR = ortho0.getMinSNR();
  double maxSNR = ortho0.getMaxSNR();

  if (myid==0) {
    std::cout << "Found minSNR=" << minSNR
	      << " maxSNR=" << maxSNR << std::endl;
  }

  if (maxSNR < minSNR )  minSNR = maxSNR * 1.0e-2;

  if (vm.count("minSNR")) {
    if (minSNR < minSNR0)  minSNR = minSNR0;
  }
  
  if (LOG) {
    minSNR = log(minSNR);
    maxSNR = log(maxSNR);
  }
  
  double dSNR = (maxSNR - minSNR)/(NSNR - 1);

  if (myid==0) {
    std::cout << "Using minSNR=" << minSNR
	      << " maxSNR=" << maxSNR
	      << " dSNR=" << dSNR << std::endl;
  }


  using OrthoCoefs = std::vector<Vector>;
  std::vector<OrthoCoefs> ac_cos(coefs.size()), ac_sin(coefs.size());
  for (int j=0; j<coefs.size(); j++) {
    ac_cos[j].resize(mmax+1), ac_sin[j].resize(mmax+1);
    for (int mm=0; mm<=mmax; mm++) {
      ac_cos[j][mm].setsize(0, norder-1);
      ac_sin[j][mm].setsize(0, norder-1);
    }
  }

  for (int nsnr=0; nsnr<NSNR; nsnr++) {

    // Assign the snr value
    //
    double snr = minSNR + dSNR*nsnr;
    if (LOG) snr = exp(snr);
    
    if (myid==0) {
      std::cout << "Computing SNR=" << snr;
      if (Hall) std::cout << " using Hall smoothing . . . " << std::endl;
      else      std::cout << " using truncation . . . "     << std::endl;
    }
    
    // Get the snr trimmed coefficients
    //
    if (myid==0) {
      std::cout << std::endl << "Trimming coefficients . . ." << std::endl;
      progress = boost::make_shared<boost::progress_display>(coefs.size());
    }

    for (int j=0; j<coefs.size(); j++) {
   
      for (int mm=0; mm<=mmax; mm++) {
	ortho0.set_coefs(mm, coefs[j]->coefC[mm], coefs[j]->coefS[mm]);
      }

      ortho0.get_trimmed(snr, ac_cos[j], ac_sin[j]);
      if (myid==0) ++(*progress);
    }
    
    // Particle loop again for KL
    //
    p = psp->GetParticle();
    int icnt = 0, ibnch = 0;

    if (myid==0) std::cout << std::endl
			   << "Computing KL for subsamples . . . "
			   << std::endl;

				// KL values, density workspace
    std::vector<double> KL(coefs.size(), 0.0), DD(coefs.size());

    unsigned good = 0, bad = 0;
    double tmas = 0.0;

    if (myid==0) {
      progress = boost::make_shared<boost::progress_display>(nbunch0*nbunch1);
    }

    do {
      if (myid==0) ++(*progress);
				// Done processing
      if (icnt >= nbunch0*nbunch1) {
	if (myid==0)
	  std::cout << "Finished KL subsamples with n=" << icnt
		    << " N=" << nbod << std::endl;
	break;
      }
				// Start a new bunch?
      if (icnt > 0 and icnt % nbunch1 == 0) ibnch++;
      
				// Particle accumulation
      if (icnt % numprocs == myid) {

				// Compute density basis for each particle
	double R   = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
	double phi = atan2(p->pos(1), p->pos(0));
	double z   = p->pos(2);

				// Get density grid interpolated entries
	std::fill(DD.begin(), DD.end(), 0.0);
	for (int mm=0; mm<=mmax; mm++) {
	  for (int nn=0; nn<norder; nn++) {
	    double dC, dS;
	    ortho0.getDensSC(mm, nn, R, z, dC, dS);
				// Sum over all subsamples
	    for (int j=0; j<coefs.size(); j++) {
	      if (j==ibnch) {
		DD[j] += coefs[j]->coefC[mm][nn]*dC*cos(phi*mm);
		if (mm) DD[j] += coefs[j]->coefS[mm][nn]*dS*sin(phi*mm);
	      } else {
		DD[j] += ac_cos[j][mm][nn]*dC*cos(phi*mm);
		if (mm) DD[j] += ac_sin[j][mm][nn]*dS*sin(phi*mm);
	      }
	    }
	  }
	}

	for (int j=0; j<coefs.size(); j++) {
	  if (j==ibnch) continue;
	  if (Ndens) {
	    if (KDdens[icnt]>0.0 and DD[j]>0.0) {
	      KL[ibnch] += p->mass() * log(KDdens[icnt]/DD[j]);
	      good++;
	      if (false) {
		static int jcnt = 0;
		if (myid==0 and j==0 and jcnt<30) {
		  std::cout << "DENS: "
			    << std::setw( 6) << p->indx()
			    << std::setw(18) << KDdens[icnt]
			    << std::setw(18) << DD[j]
			    << std::endl;
		  jcnt++;
		}
	      }
	    } else {
	      bad++;
	    }
	  } else {
	    if (DD[ibnch]>0.0 and DD[j]>0.0) {
	      KL[ibnch] += p->mass() * log(DD[ibnch]/DD[j]);
	      good++;
	    } else {
	      bad++;
	    }
	  }
	}

	tmas += p->mass();
      }
      // END: parallelized work
      
      p = psp->NextParticle();
      icnt++;
    } while (p);
    
    // For diagnostic output
    //
    if (myid) {
      MPI_Reduce(&good, 0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&bad , 0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&tmas, 0, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(MPI_IN_PLACE, &good, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &bad , 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &tmas, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Sum reduce KL from all processes
    //
    if (myid) {
      MPI_Reduce(KL.data(), 0, coefs.size(),
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Root node records results
    //
    if (myid==0) {
      double ratio = static_cast<double>(bad)/good;
      std::cout << std::endl << "Bad/good density counts ["
		<< bad << "/" << good << "=" << ratio << "]" << std::endl;
      double corr = log(1.0 + ratio);

      MPI_Reduce(MPI_IN_PLACE, KL.data(), coefs.size(),
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      out << std::setw(18) << snr << std::setw(18)
	  << std::accumulate(KL.begin(), KL.end(), 0.0) * ampfac/tmas + corr
	  << std::setw(18) << ratio
	  << std::setw(18) << corr
	  << std::endl;
    }
    // END: root node records results

  }
  // END: SNR loop
      
  // Blank line between stanzas
  //
  if (myid==0) out << std::endl;


  // Clean up and exit
  //
  MPI_Finalize();

  return 0;
}

