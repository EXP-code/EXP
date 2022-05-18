/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  k-d density test for cylinder
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
 *  Output matrix file for Gnuplot
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

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <EmpCylSL.H>
#include <foarray.H>
#include <KDtree.H>

#include <localmpi.H>
#include <Progress.H>
#include <cxxopts.H>
#include <libvars.H>

#include <yaml-cpp/yaml.h>	// YAML support


int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, rscale, minSNR0, Hexp, ROUT, ZOUT;
  int NICE, LMAX, NMAX, NSNR, indx, nbunch, Ndens, NOUT, NPHI;
  std::string cachefile, dir("./"), cname, prefix, fileType, psfile;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "K-D density estimate test" << std::endl
       << std::string(60, '-') << std::endl << std::endl
       << "Allowed options";
  
  cxxopts::Options options("kdtest", sout.str());

  options.add_options()
    ("h,help", "Print this help message")
    ("K,Ndens", "KD density estimate count (use 0 for expansion estimate)",
     cxxopts::value<int>(Ndens)->default_value("32"))
    ("NICE", "system priority",
     cxxopts::value<int>(NICE)->default_value("0"))
    ("NOUT", "Number of grid points for surface output",
     cxxopts::value<int>(NOUT)->default_value("40"))
    ("ROUT", "Maximum radius for output",
     cxxopts::value<double>(ROUT)->default_value("0.05"))
    ("ZOUT", "Maximum height for output",
     cxxopts::value<double>(ZOUT)->default_value("0.01"))
    ("NPHI", "Number of azimuthal bins",
     cxxopts::value<int>(NPHI)->default_value("16"))
    ("prefix", "Filename prefix",
     cxxopts::value<std::string>(prefix)->default_value("kdtest"))
    ("runtag", "Phase space file",
     cxxopts::value<std::string>(runtag)->default_value("run1"))
    ("outdir", "Output directory path",
     cxxopts::value<std::string>(outdir)->default_value("."))
    ("psfile", "phase-space file name name",
     cxxopts::value<std::string>(psfile))
    ("cachefile", "cachefile name",
     cxxopts::value<std::string>(cachefile)->default_value(".eof.cache.file"))
    ("cname","component name",
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

    std::ifstream in(cachefile);
    if (!in) {
      std::cerr << "Error opening cachefile named <" 
		<< cachefile << "> . . ."
		<< std::endl
		<< "I will build <" << cachefile
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

	std::cout << "idens=" << idens << std::endl;

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

				// Create expansion
				//
  EmpCylSL ortho(NMAX, LMAX, mmax, norder, rscale, vscale, nodd, cachefile);
    
  if (ortho.read_cache()==0) {
    std::cout << "Could not read cache file <" << cachefile << ">"
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
  // Phase space
  // ==================================================

#ifdef DEBUG
  std::cout << "Begin phase -space loop" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif	      

  // ==================================================
  // Open output file
  // ==================================================

  std::ofstream out;
  out.open(prefix + ".out");
  if (!out) {
    std::cerr << "Error opening output file <" << prefix + ".out" << ">" << std::endl;
    exit(-2);
  }
  
  // ==================================================
  // Open PSP file
  // ==================================================

  std::vector<std::string> file1 = {psfile};

  PR::PRptr reader = PR::ParticleReader::createReader
    (fileType, file1, myid, true);
  
  double tnow = reader->CurrentTime();
  if (myid==0) std::cout << "Beginning partition [time=" << tnow
			 << ", index=" << indx << "] . . . "  << flush;
  
  reader->SelectType(cname);

  int nbod = reader->CurrentNumber();

  std::shared_ptr<progress::progress_display> progress;
  if (myid==0) {
    std::cout << std::endl
	      << "Accumulating particle positions . . . "
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(nbod);
  }

  auto p = reader->firstParticle();
  int icnt = 0;

  ortho.setup_accumulation();

  do {
    if (myid==0) ++(*progress);

    if (icnt++ % numprocs == myid) {
      double R   = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
      double phi = atan2(p->pos[1], p->pos[0]);
      ortho.accumulate(R, p->pos[2], phi, p->mass, p->indx, 0);
    }
    p = reader->nextParticle();
  } while (p);
  

  if (myid==0) cout << "Making coefficients for total . . . " << flush;
  ortho.make_coefficients(true);
  if (myid==0) std::cout << "done" << endl;

  // This is the kd- NN density estimate; skipped by default for Ndens=0
  //
  if (Ndens<=0) Ndens = 2;
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
    
  // Make the k-d tree
  //
  tree3 tree(points.begin(), points.end());
    
  // Field storage for parallel reduction
  //
  std::vector<double> kdens(NOUT*NOUT), odens(NOUT*NOUT), opot(NOUT*NOUT);
  std::fill(kdens.begin(), kdens.end(), 0.0);
  std::fill(odens.begin(), odens.end(), 0.0);
  std::fill(opot.begin(),  opot.end(),  0.0);

  if (myid==0) {
    std::cout << std::endl
	      << "Evaluating field quantities . . . "
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(NOUT*NOUT);
  }

  // Evaluate the fields
  //
  double dR = 2.0*ROUT/(NOUT-1);
  double dZ = 2.0*ZOUT/(NOUT-1);
  icnt = 0;

  for (int j=0; j<NOUT; j++) {
    double Z = -ZOUT + dZ*j;

    for (int i=0; i<NOUT; i++) {
      double R = -ROUT + dR*i;
      
      if (myid==0) ++(*progress);

      if (icnt++ % numprocs == myid) {

	double dphi = 2.0*M_PI/NPHI;
	for (int nphi=0; nphi<NPHI; nphi++) {
	  double phi = dphi*nphi;
	  auto ret   = tree.nearestN({R*cos(phi), R*sin(phi), Z}, Ndens);
	  double volume = 4.0*M_PI/3.0*std::pow(std::get<2>(ret), 3.0);
	  if (volume>0.0) kdens[j*NOUT + i] += std::get<1>(ret)/volume/NPHI;
	}

	double d, p;

	odens[j*NOUT + i] = ortho.accumulated_dens_eval(fabs(R), Z, atan2(0.0, R), d);
	ortho.accumulated_eval(fabs(R), Z, atan2(0.0, R), d, p, d, d, d);
	opot[j*NOUT + i] = p;
      }
    }
  }

  // Reduction and output
  //
  if (myid==0) {
    
    MPI_Reduce(MPI_IN_PLACE, kdens.data(), kdens.size(), MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, odens.data(), odens.size(), MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, opot.data(),  opot.size(),  MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);

    for (int j=0; j<NOUT; j++) {
      double Z = -ZOUT + dZ*j;
      
      for (int i=0; i<NOUT; i++) {
	double R = -ROUT + dR*i;
      
	out << std::setw(18) << R
	    << std::setw(18) << Z
	    << std::setw(18) << kdens[NOUT*j + i]
	    << std::setw(18) << odens[NOUT*j + i]
	    << std::setw(18) << opot [NOUT*j + i]
	    << std::endl;
      }
      out << std::endl;
    }
    std::cout << "Done" << std::endl;
  } else {
    MPI_Reduce(kdens.data(), 0, kdens.size(), MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(odens.data(), 0, odens.size(), MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(opot.data(),  0, opot.size(),  MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}

