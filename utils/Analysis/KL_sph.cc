/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Kullback-Liebler analysis for sphere
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
#include <SphereSL.H>
#include <foarray.H>

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
  int lmax, nmax;

public:

  Matrix coefs;

  CoefStruct(int lmax, int nmax) : lmax(lmax), nmax(nmax)
  {
    Matrix ret(0, lmax*(lmax+2), 1, nmax);
    ret.zero();
  }

  void sync()
  {
    for (int l=0; l<(lmax+1)*(lmax+1); l++) {
      MPI_Allreduce(MPI_IN_PLACE, &coefs[l][1], nmax,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
  int NICE, LMAX, NMAX, NSNR, indx, nbunch;
  std::string modelf, dir("./"), cname, prefix;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Kullback-Liebler analysis for spherical models" << std::endl
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
    ("NICE",                po::value<int>(&NICE)->default_value(0),
     "system priority")
    ("LMAX",                po::value<int>(&LMAX)->default_value(8),
     "Maximum harmonic order for spherical expansion")
    ("NMAX",                po::value<int>(&NMAX)->default_value(18),
     "Maximum harmonic order for spherical expansion")
    ("NSNR, N",             po::value<int>(&NSNR)->default_value(20),
     "Number of SNR evaluations")
    ("minSNR",              po::value<double>(&minSNR0),
     "minimum SNR value for loop output")
    ("rscale",              po::value<double>(&rscale)->default_value(0.067),
     "Radial scale for coordinate mapping (cmap)")
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
    ("modelfile",
     po::value<std::string>(&modelf)->default_value("SLGridSph.model"),
     "halo model file name")
    ("cname",
     po::value<std::string>(&cname)->default_value("dark halo"),
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
  // PSP input stream
  // ==================================================

  int iok = 1;
  std::ostringstream s1;
  
  s1.str("");		// Clear stringstream
  if (SPL) s1 << "SPL.";
  else     s1 << "OUT.";
  s1 << runtag << "."<< std::setw(5) << std::setfill('0') << indx;
      
				// Check for existence of next file
  std:string file = dir + s1.str();
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

  PSPptr psp;

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
      
  // ==================================================
  // Make SL expansion
  // ==================================================

  SphericalModelTable halo(modelf);
  SphereSL::mpi  = true;
  SphereSL::NUMR = 4000;
  SphereSL::HEXP = Hexp;

  int nbod = psp->GetNamed(cname)->comp.nbod;
  int nprt = std::floor(sqrt(nbod));

  SphereSL ortho0(&halo, LMAX, NMAX, 1, rscale, true, nprt);
  SphereSL ortho1(&halo, LMAX, NMAX, 1, rscale);

  if (myid==0) std::cout << std::endl
			 << "Accumulating particle positions . . . "
			 << std::endl;

  // Zero out coefficients to prepare for a new expansion
  //
  ortho0.reset_coefs();

  boost::shared_ptr<boost::progress_display> progress;
  if (myid==0) {
    progress = boost::make_shared<boost::progress_display>(nbod);
  }

  SParticle *p = psp->GetParticle();
  int icnt = 0;
  do {
    if (myid==0) ++(*progress);

    if (icnt++ % numprocs == myid) {
      ortho0.accumulate(p->pos(0), p->pos(1), p->pos(2), p->mass());
    }
    p = psp->NextParticle();
  } while (p);
  
    
  //------------------------------------------------------------ 
      
  if (myid==0) std::cout << std::endl
			 << "Making coefficients for total . . . "
			 << std::flush;
  ortho0.make_coefs();
  ortho0.make_covar();

  if (myid==0) std::cout << "done" << endl;
  
  if (myid==0) std::cout << std::endl
			 << "Accumulating particle positions for subsamples . . . "
			 << std::endl;

				// Size of bunch
  int nbunch1 = std::floor(sqrt(nbod));
  if (nbunch>0) nbunch1 = nbod/nbunch;
				// Number of bunches
  int nbunch0 = nbod/nbunch1;

  p = psp->GetParticle();
  icnt = 0;
    
  std::vector<boost::shared_ptr<CoefStruct>> coefs;

  if (myid==0) {
    progress = boost::make_shared<boost::progress_display>(nbunch0*nbunch1);
  }

  do {
    if (myid==0) ++(*progress);
				// Done processing
    if (icnt >= nbunch0*nbunch1) {
      if (myid==0)
	std::cout << "Finished processing subsamples with " << icnt
		  << "/" << nbod << std::endl;
      break;
    }
      
				// Start a new bunch?
    if (icnt % nbunch1 == 0) {
      if (coefs.size()) {
	ortho1.make_coefs();
	coefs.back()->coefs = ortho1.retrieve_coefs();
	coefs.back()->sync();
      }
      coefs.push_back(boost::make_shared<CoefStruct>(LMAX, NMAX));
      ortho1.reset_coefs();
    }
    
				// Particle accumulation
    if (icnt++ % numprocs == myid) {
      ortho1.accumulate(p->pos(0), p->pos(1), p->pos(2), p->mass());
    }
    p = psp->NextParticle();
  } while (p);
  
  //------------------------------------------------------------ 
      
  if (myid==0) cout << std::endl
		    << "Making coefficients for bunch ["
		    << coefs.size() << "/" << nbunch0 << "] " << flush;

  ortho1.make_coefs();
  coefs.back()->coefs = ortho1.retrieve_coefs();
  coefs.back()->sync();

  if (myid==0) std::cout << "done" << endl;
  
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

  for (int nsnr=0; nsnr<NSNR; nsnr++) {

    // Assign the snr value
    //
    double snr = minSNR + dSNR*nsnr;
    if (LOG) snr = exp(snr);
    
    if (myid==0) {
      std::cout << "Computing SNR=" << snr;
      if (Hall) std::cout << " using Hall smoothing." << std::endl;
      else      std::cout << " using truncation." << std::endl;
    }
    
    // Get the snr trimmed coefficients
    //
    std::vector<Matrix> coefs1(coefs.size());

    if (myid==0) {
      std::cout << std::endl << "Trimming coefficients . . ." << std::endl;
      progress = boost::make_shared<boost::progress_display>(coefs.size());
    }

    for (int j=0; j<coefs.size(); j++) {
      ortho0.install_coefs(coefs[j]->coefs);
      coefs1[j] = ortho0.get_trimmed(snr, Hall);
      if (myid==0) ++(*progress);
    }
    
    // Particle loop again for KL
    //
    p = psp->GetParticle();
    int icnt = 0, ibnch = 0;

    if (myid==0) std::cout << std::endl
			   << "Computing KL for subsamples . . ."
			   << std::endl;

				// KL values, density workspace
    std::vector<double> KL(coefs.size(), 0.0), DD(coefs.size());

    unsigned good = 0, bad = 0;

    if (myid==0) {
      progress = boost::make_shared<boost::progress_display>(nbunch0*nbunch1);
    }

    do {
      if (myid==0) ++(*progress);
				// Done processing
      if (icnt >= nbunch0*nbunch1) {
	if (myid==0)
	  std::cout << "Finished KL subsamples with " << icnt
		    << "/" << nbod << std::endl;
	break;
      }
				// Start a new bunch?
      if (icnt % nbunch1 == 0) ibnch++;
      
				// Particle accumulation
      if (icnt++ % numprocs == myid) {

				// Compute density for each particle
	double r     = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1) + p->pos(2)*p->pos(2));
	double costh = p->pos(2)/(r + 1.0e-18);
	double phi   = atan2(p->pos(1), p->pos(1));

	for (int j=0; j<coefs.size(); j++) {
	  double t0, t1, t2, t3;
	  ortho1.install_coefs(coefs1[j]);
	  ortho1.dens_pot_eval(r, costh, phi, t0, t1, t2, t3);
	  DD[j] = t0 + t1;
	}

	for (int j=0; j<coefs.size(); j++) {
	  if (DD[ibnch]>0.0 and DD[j]>0.0) {
	    KL[j] += p->mass() * log(DD[ibnch]/DD[j]);
	    good++;
	  } else {
	    bad++;
	  }
	}

      }
      p = psp->NextParticle();
    } while (p);
    
    // For diagnostic output
    //
    if (myid) {
      MPI_Reduce(&good, 0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&bad , 0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(MPI_IN_PLACE, &good, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &bad , 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Sum reduce KL from all processes
    //
    if (myid) 
      MPI_Reduce(KL.data(), 0, coefs.size(),
	      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0) {
      std::cout << std::endl << "Bad/good density counts ["
		<< bad << "/" << good << "="
		<< static_cast<double>(bad)/good << "]" << std::endl;

      MPI_Reduce(MPI_IN_PLACE, KL.data(), coefs.size(),
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      out << std::setw(18) << snr << std::setw(18)
	  << std::accumulate(KL.begin(), KL.end(), 0.0)
	  << std::endl;
    }

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

