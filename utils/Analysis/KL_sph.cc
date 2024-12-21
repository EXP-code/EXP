/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Kullback-Leibler analysis for sphere
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
				// Progress barp
#include <Progress.H>
                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <SphSL.H>
#include <KDtree.H>
#include <foarray.H>
#include <cxxopts.H>		// Command-line parsing
				// Library support
#include <localmpi.H>
#include <libvars.H>

using namespace __EXP__;	// Reference to n-body globals

#include <yaml-cpp/yaml.h>	// YAML support

  
// Helper class
//
class CoefStruct
{
private:
  int lmax, nmax;

public:

  Eigen::MatrixXd coefs;

  CoefStruct(int lmax, int nmax) : lmax(lmax), nmax(nmax)
  {
    Eigen::MatrixXd ret((lmax+1)*(lmax+1), nmax);
    ret.setZero();
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
    coefs /= norm;
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
  std::string modelf, dir("./"), cname, fileType, psfile, prefix;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Kullback-Leibler analysis for spherical models" << std::endl
       << std::string(60, '-') << std::endl << std::endl;
  
  cxxopts::Options options(argv[0], sout.str());

  options.add_options()
    ("h,help", "Print this help message")
    ("v,verbose", "Verbose and diagnostic output for covariance computation")
    ("t,truncate", "Use Truncate method for SNR trimming rather than the default Hall")
    ("d,debug", "Debug max values")
    ("LOG", "log scaling for SNR")
    ("Hall", "use Hall smoothing for SNR trim")
    ("K,Ndens", "KD density estimate count (use 0 for expansion estimate)",
     cxxopts::value<int>(Ndens)->default_value("32"))
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("psfile", "Phase-space file",
     cxxopts::value<std::string>(psfile))
    ("NICE", "system priority",
     cxxopts::value<int>(NICE)->default_value("0"))
    ("LMAX", "Maximum harmonic order for spherical expansion",
     cxxopts::value<int>(LMAX)->default_value("8"))
    ("NMAX", "Maximum harmonic order for spherical expansion",
     cxxopts::value<int>(NMAX)->default_value("18"))
    ("N,NSNR", "Number of SNR evaluations",
     cxxopts::value<int>(NSNR)->default_value("20"))
    ("minSNR", "minimum SNR value for loop output",
     cxxopts::value<double>(minSNR0))
    ("rscale", "Radial scale for coordinate mapping (cmap)",
     cxxopts::value<double>(rscale)->default_value("0.067"))
    ("Hexp", "default Hall smoothing exponent",
     cxxopts::value<double>(Hexp)->default_value("1.0"))
    ("prefix", "Filename prefix",
     cxxopts::value<string>(prefix)->default_value("KLsph"))
    ("runtag", "Phase space file",
     cxxopts::value<string>(runtag)->default_value("run1"))
    ("outdir", "Output directory path",
     cxxopts::value<string>(outdir)->default_value("."))
    ("indx", "PSP index",
     cxxopts::value<int>(indx)->default_value("0"))
    ("nbunch", "Desired bunch size (default: sqrt(nbod) if value is < 0)",
     cxxopts::value<int>(nbunch)->default_value("-1"))
    ("modelfile", "halo model file name",
     cxxopts::value<std::string>(modelf)->default_value("SLGridSph.model"))
    ("cname", "component name",
     cxxopts::value<std::string>(cname)->default_value("dark halo"))
    ;
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    if (myid==0) std::cout << std::endl << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }

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

  PR::PRptr reader = PR::ParticleReader::createReader
    (fileType, {psfile}, myid, true);
  
  double tnow = reader->CurrentTime();
  if (myid==0) std::cout << "Beginning partition [time=" << tnow
			 << ", index=" << indx << "] . . . "  << flush;
  
  
  reader->SelectType(cname);
      
  // ==================================================
  // Make SL expansion
  // ==================================================

  auto halo = std::make_shared<SphericalModelTable>(modelf);

  SphSL::mpi  = true;
  SphSL::NUMR = 4000;
  SphSL::HEXP = Hexp;

  int nbod = reader->CurrentNumber();
  int nprt = std::floor(sqrt(nbod));

  SphSL ortho0(halo, LMAX, NMAX, 1, rscale, true, nprt);
  SphSL ortho1(halo, LMAX, NMAX, 1, rscale);

  if (myid==0) std::cout << std::endl
			 << "Accumulating particle positions . . . "
			 << std::endl;

  // Zero out coefficients to prepare for a new expansion
  //
  ortho0.reset_coefs();

  std::vector<double> KDdens;

  std::shared_ptr<progress::progress_display> progress;
  if (myid==0) {
    progress = std::make_shared<progress::progress_display>(nbod);
  }

  auto p = reader->firstParticle();
  int icnt = 0;
  do {
    if (myid==0) ++(*progress);

    if (icnt++ % numprocs == myid) {
      ortho0.accumulate(p->pos[0], p->pos[1], p->pos[2], p->mass);
    }
    p = reader->nextParticle();
  } while (p);
  
    
  // This is the kd- NN density estimate; skipped by default for Ndens=0
  //
  if (Ndens) {
    if (myid==0) std::cout << "Computing KD density estimate for " << nbod
			   << " points" << std::endl;

    using point3 = KDtree::point <double, 3>;
    using tree3  = KDtree::kdtree<double, 3>;

    std::vector<point3> points;

    // Every node needs to make the tree (a parallel share could be
    // implemented in KDtree.H)
    //
    double KDmass = 0.0;
    for (auto part=reader->firstParticle(); part!=0; part=reader->nextParticle()) {
      KDmass += part->mass;
      points.push_back(point3({part->pos[0], part->pos[1], part->pos[2]}, part->mass));
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
      
  if (myid==0) std::cout << std::endl
			 << "Making coefficients for total . . . "
			 << std::flush;
  ortho0.make_coefs();
  ortho0.make_covar(verbose);

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

				// Reset the particle iterator
  p = reader->firstParticle();
  icnt = 0;
    
  std::vector<std::shared_ptr<CoefStruct>> coefs;

  if (myid==0) {
    progress = std::make_shared<progress::progress_display>(nbunch0*nbunch1);
  }

  double curMass = 0.0;

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
	coefs.back()->sync(curMass);
      }
      coefs.push_back(std::make_shared<CoefStruct>(LMAX, NMAX));
      ortho1.reset_coefs();
      curMass = 0.0;
    }
    
				// Particle accumulation is spread
				// between nodes
    if (icnt++ % numprocs == myid) {
      ortho1.accumulate(p->pos[0], p->pos[1], p->pos[2], p->mass);
      curMass += p->mass;
    }
    p = reader->nextParticle();
  } while (p);
  
  //------------------------------------------------------------ 
      
  if (myid==0) cout << std::endl
		    << "Making coefficients for bunch ["
		    << coefs.size() << "/" << nbunch0 << "] " << flush;

  ortho1.make_coefs();
  coefs.back()->coefs = ortho1.retrieve_coefs();
  coefs.back()->sync(curMass);

  if (myid==0) std::cout << "done" << endl;
  
  // This is a debug test
  if (myid==0) {
    std::ofstream test(prefix + ".coeftest");
    for (int n=0; n<NMAX; n++) {
      test << std::setw(4) << n;
      for (auto & c : coefs)
	test << std::setw(18) << c->coefs(0, n);
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
    std::vector<Eigen::MatrixXd> coefs1(coefs.size());

    if (myid==0) {
      std::cout << std::endl << "Trimming coefficients . . ." << std::endl;
      progress = std::make_shared<progress::progress_display>(coefs.size());
    }

    for (int j=0; j<coefs.size(); j++) {
      ortho0.install_coefs(coefs[j]->coefs);
      coefs1[j] = ortho0.get_trimmed(snr, 1.0, Hall);
      if (myid==0) ++(*progress);
    }
    
    // Particle loop again for KL
    //
    p = reader->firstParticle();
    int icnt = 0, ibnch = 0;

    if (myid==0) std::cout << std::endl
			   << "Computing KL for subsamples . . ."
			   << std::endl;

				// KL values, density workspace
    std::vector<double> KL(coefs.size(), 0.0), DD(coefs.size());

    unsigned good = 0, bad = 0;
    double tmas = 0.0;

    if (myid==0) {
      progress = std::make_shared<progress::progress_display>(nbunch0*nbunch1);
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
      if (icnt > 0 and icnt % nbunch1 == 0) ibnch++;
      
				// Particle accumulation
      if (icnt % numprocs == myid) {

				// Compute density for each particle
	double r     = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1] + p->pos[2]*p->pos[2]);
	double costh = p->pos[2]/(r + 1.0e-18);
	double phi   = atan2(p->pos[1], p->pos[1]);

				// Make the density for each bunch
	for (int j=0; j<coefs.size(); j++) {
	  double t0, t1, t2, t3;
	  if (ibnch == j)
	    ortho1.install_coefs(coefs[j]->coefs);
	  else
	    ortho1.install_coefs(coefs1[j]);

	  ortho1.dens_pot_eval(r, costh, phi, t0, t1, t2, t3);
	  DD[j] = t0 + t1;
	}

	for (int j=0; j<coefs.size(); j++) {
	  if (j==ibnch) continue;
	  if (Ndens) {
	    if (KDdens[icnt]>0.0 and DD[j]>0.0) {
	      KL[ibnch] += p->mass * log(KDdens[icnt]/DD[j]);
	      good++;
	      if (false) {
		static int jcnt=0;
		if (myid==0 and j==0) {
		  if (jcnt++ < 30) {
		    std::cout << "DENS: "
			      << std::setw(18) << KDdens[icnt]
			      << std::setw(18) << DD[j]
			      << std::endl;
		  }
		}
	      }
	    } else {
	      bad++;
	    }
	  } else {
	    if (DD[ibnch]>0.0 and DD[j]>0.0) {
	      KL[ibnch] += p->mass * log(DD[ibnch]/DD[j]);
	      good++;
	    } else {
	      bad++;
	    }
	  }
	}

	tmas += p->mass;
      }
      // END: parallelized work
      
      p = reader->nextParticle();
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
    if (myid) 
      MPI_Reduce(KL.data(), 0, coefs.size(),
	      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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

