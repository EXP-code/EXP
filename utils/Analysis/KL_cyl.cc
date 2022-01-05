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
#include <random>
#include <queue>
#include <map>

				// Eigen3
#include <Eigen/Eigen>

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <EmpCylSL.H>
#include <foarray.H>
#include <KDtree.H>
#include <Progress.H>
#include <cxxopts.H>

#include <global.H>
#include <localmpi.H>

#include <yaml-cpp/yaml.h>	// YAML support

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
  std::string CACHEFILE, dir("./"), cname, prefix, fileType, filePrefix;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Kullback-Leibler analysis for cylindrical models" << std::endl
       << std::string(60, '-') << std::endl << std::endl;
  
  cxxopts::Options options(argv[0], sout.str());

  options.add_options()
    ("H,help", "Print this help message")
    ("v,verbose", "Verbose and diagnostic output for covariance computation")
    ("t,truncate", "Use Truncate method for SNR trimming rather than the default Hall")
    ("debug", "Debug max values")
    ("LOG", "log scaling for SNR")
    ("Hall", "use Hall smoothing for SNR trim")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("P,prefix", "prefix for phase-space files",
     cxxopts::value<std::string>(filePrefix)->default_value("OUT"))
    ("K,Ndens", "KD density estimate count (use 0 for expansion estimate)",
     cxxopts::value<int>(Ndens)->default_value("32"))
    ("NICE", "system priority",
     cxxopts::value<int>(NICE)->default_value("0"))
    ("LMAX", "Maximum harmonic order for spherical expansion",
     cxxopts::value<int>(LMAX)->default_value("36"))
    ("NSNR, N", "Number of SNR evaluations",
     cxxopts::value<int>(NSNR)->default_value("20"))
    ("minSNR", "minimum SNR value for loop output",
     cxxopts::value<double>(minSNR0))
    ("Hexp", "default Hall smoothing exponent",
     cxxopts::value<double>(Hexp)->default_value("1.0"))
    ("prefix", "Filename prefix",
     cxxopts::value<string>(prefix)->default_value("crossval"))
    ("runtag", "Phase space file",
     cxxopts::value<string>(runtag)->default_value("run1"))
    ("outdir", "Output directory path",
     cxxopts::value<string>(outdir)->default_value("."))
    ("indx", "PSP index",
     cxxopts::value<int>(indx)->default_value("0"))
    ("nbunch", "Desired bunch size (default: sqrt(nbod) if value is < 0)",
     cxxopts::value<int>(nbunch)->default_value("-1"))
    ("dir,d", "directory for SPL files",
     cxxopts::value<std::string>(dir))
    ("ignore", "rebuild EOF grid if input parameters do not match the cachefile",
     cxxopts::value<bool>(ignore)->default_value("false"))
    ("cachefile", "cachefile name",
     cxxopts::value<std::string>(CACHEFILE)->default_value(".eof.cache.file"))
    ("cname", "component name",
     cxxopts::value<std::string>(cname)->default_value("star disk"))
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
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  int mmax, numx, numy, norder, cmapr, cmapz, nodd=-1;
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
	NMAX    = node["nmax"  ].as<int>();
	norder  = node["norder"].as<int>();
	DENS    = node["dens"  ].as<bool>();
	if (node["nodd"])
	  nodd  = node["nodd"  ].as<int>();
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
  EmpCylSL::PCAVAR      = true;
  EmpCylSL::PCADRY      = true;

				// Create expansion
				//
  EmpCylSL ortho0(NMAX, LMAX, mmax, norder, rscale, vscale, nodd, CACHEFILE);
  EmpCylSL ortho1(NMAX, LMAX, mmax, norder, rscale, vscale, nodd, CACHEFILE);
    
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
  auto file1 = PR::ParticleReader::fileNameCreator
    (fileType, indx, myid, dir, runtag);

  std::ifstream in(file);
  if (!in) {
    if (myid==0) 
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

  PR::PRptr reader = PR::ParticleReader::createReader
    (fileType, file1, myid, true);

  double tnow = reader->CurrentTime();
  if (myid==0) std::cout << "Beginning partition [time=" << tnow
			 << ", index=" << indx << "] . . . "  << flush;
  
  reader->SelectType(cname);
      
  int nbod = reader->CurrentNumber();

  std::vector<double> KDdens;

  std::shared_ptr<progress::progress_display> progress;
  if (myid==0) {
    std::cout << std::endl
	      << "Accumulating particle positions . . . "
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(nbod);
  }

  ortho0.setup_accumulation();
  ortho0.setHall("test", nbod);

  auto p = reader->firstParticle();
  int icnt = 0;
  do {
    if (myid==0) ++(*progress);

    if (icnt++ % numprocs == myid) {
      double R   = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
      double phi = atan2(p->pos[1], p->pos[0]);
      ortho0.accumulate(R, p->pos[2], phi, p->mass, p->indx, 0, 0, true);
      //                                                     ^  ^  ^
      //                                                     |  |  |
      // Thread id ------------------------------------------+  |  |
      // Level -------------------------------------------------+  |
      // Compute covariance ---------------------------------------+
    }
    p = reader->nextParticle();
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
      
  if (myid==0) cout << "Making coefficients for total . . . " << flush;
  ortho0.make_coefficients(true);
  ortho0.pca_hall(true, false);
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
      coefs.push_back(std::make_shared<CoefStruct>(mmax, norder));
      ortho1.setup_accumulation();
      curMass = 0.0;
    }
    
				// Particle accumulation is spread
				// between nodes
    if (icnt++ % numprocs == myid) {
      
      double R   = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
      double phi = atan2(p->pos[1], p->pos[0]);
      ortho1.accumulate(R, p->pos[2], phi, p->mass, p->indx, 0, 0, false);
      //                                                     ^  ^  ^
      //                                                     |  |  |
      // Thread id ------------------------------------------+  |  |
      // Level -------------------------------------------------+  |
      // Compute covariance ---------------------------------------+

      curMass += p->mass;	// Accumulate mass per bunch
    }
    p = reader->nextParticle();
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


  using OrthoCoefs = std::vector<Eigen::VectorXd>;
  std::vector<OrthoCoefs> ac_cos(coefs.size()), ac_sin(coefs.size());
  for (int j=0; j<coefs.size(); j++) {
    ac_cos[j].resize(mmax+1), ac_sin[j].resize(mmax+1);
    for (int mm=0; mm<=mmax; mm++) {
      ac_cos[j][mm].resize(norder);
      ac_sin[j][mm].resize(norder);
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
      progress = std::make_shared<progress::progress_display>(coefs.size());
    }

    for (int j=0; j<coefs.size(); j++) {
   
      for (int mm=0; mm<=mmax; mm++) {
	if (mm==0)
	  ortho0.set_coefs(mm, coefs[j]->coefC[mm], coefs[j]->coefS[mm], true);
	else
	  ortho0.set_coefs(mm, coefs[j]->coefC[mm], coefs[j]->coefS[mm], false);
      }

      ortho0.get_trimmed(snr, ac_cos[j], ac_sin[j]);
      if (myid==0) ++(*progress);
    }
    
    // Reset particle loop again for KL
    //
    p = reader->firstParticle();
    int icnt = 0, ibnch = 0;

    // KL values, density workspace
    //
    std::vector<double> KL(coefs.size(), 0.0), DD(coefs.size());

    unsigned good = 0, bad = 0;
    double tmas = 0.0;

    if (myid==0) {
      std::cout << std::endl
		<< "Computing KL for subsamples . . . "
		<< std::endl;
      progress = std::make_shared<progress::progress_display>(nbunch0*nbunch1);
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
	double R   = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	double phi = atan2(p->pos[1], p->pos[0]);
	double z   = p->pos[2];
	
				// Get density grid interpolated entries
	std::fill(DD.begin(), DD.end(), 0.0);
	for (int mm=0; mm<=mmax; mm++) {
	  for (int nn=0; nn<norder; nn++) {
	    double dC, dS;
	    ortho0.getDensSC(mm, nn, R, z, dC, dS);
				// Sum over all subsamples
	    for (int j=0; j<coefs.size(); j++) {
	      DD[j] += ac_cos[j][mm][nn]*dC*cos(phi*mm);
	      if (mm) DD[j] += ac_sin[j][mm][nn]*dS*sin(phi*mm);
	    }
	  }
	}

	for (int j=0; j<coefs.size(); j++) {
	  if (j==ibnch) continue;
	  if (Ndens) {
	    if (KDdens[icnt]>0.0 and DD[j]>0.0) {
	      KL[ibnch] += p->mass * log(KDdens[icnt]/DD[j]);
	      good++;
	      if (false) {
		static int jcnt = 0;
		if (myid==0 and j==0 and jcnt<100) {
		  std::cout << "DENS: "
			    << std::setw( 6) << p->indx
			    << std::setw(18) << KDdens[icnt]
			    << std::setw(18) << DD[0]
			    << std::endl;
		  jcnt++;
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

