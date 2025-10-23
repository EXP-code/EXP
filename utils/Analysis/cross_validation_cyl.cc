/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Cross validation analysis for cylinder
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
 *  MDW 11/28/08
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
#include <memory>
#include <queue>
#include <map>
				// Eigen3
#include <Eigen/Eigen>

#include "Progress.H"

                                // System libs
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include "numerical.H"
#include "ParticleReader.H"
#include "interp.H"
#include "massmodel.H"
#include "EmpCylSL.H"
#include "foarray.H"

#include "libvars.H"
#include "cxxopts.H"
#include "localmpi.H"

#include <yaml-cpp/yaml.h>	// YAML support
  
// Globals
//
extern double plgndr(int ll, int mm, double x);

double Ylm_fac(int ll, int mm)
{
  mm = abs(mm);
  return sqrt( (2.0*ll+1)/(4.0*M_PI) ) *
    exp(0.5*(lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm)));
}

double Xi(double R, double Rp, double z, double zp)
{
  const double tol = 1.0e-18;
  if (R < 1.0e-18 or Rp < 1.0e-18)
    return std::numeric_limits<double>::infinity();

  double q = R/Rp;
  if (q < 1.0e-18 or 1.0/q < 1.0e-18)
    return std::numeric_limits<double>::infinity();

  return 0.5*(q + 1.0/q + (z - zp)*(z - zp)/(R*Rp));
}

// Fortran binding for the half-integral Legendre functions of the
// second kind
//
extern "C" {
  int dtorh1_(double* z, int* m, int* nmax,
	      double* pl, double* ql, int* newn);
}
  

int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, RMAX, rscale, minSNR;
  int NICE, LMAX, NMAX, NSNR, NPART, init, num;
  std::string CACHEFILE, modelf, cname, psfiles, delim;
  std::string prefix, table_cache, fileType;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Cross-validation analysis for cylindrical models" << std::endl
       << std::string(60, '-') << std::endl;
  
  cxxopts::Options options(argv[0], sout.str());

  options.add_options()
    ("h,help", "Print this help message")
    ("v,verbose", "Verbose and diagnostic output for covariance computation")
    ("LOG", "log scaling for SNR")
    ("Hall", "use Hall smoothing for SNR trim")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("NICE", "system priority",
     cxxopts::value<int>(NICE)->default_value("0"))
    ("RMIN", "minimum radius for output",
     cxxopts::value<double>(RMIN)->default_value("0.0"))
    ("RSCALE", "coordinate mapping scale factor",
     cxxopts::value<double>(rscale)->default_value("0.067"))
    ("RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("2.0"))
    ("LMAX", "Maximum harmonic order for spherical expansion",
     cxxopts::value<int>(LMAX)->default_value("4"))
    ("NMAX", "Maximum radial order for spherical expansion",
     cxxopts::value<int>(NMAX)->default_value("12"))
    ("NPART", "Jackknife partition number for testing (0 means off, use standard eval)",
     cxxopts::value<int>(NPART)->default_value("0"))
    (" N,NSNR", "Number of SNR evaluations",
     cxxopts::value<int>(NSNR)->default_value("20"))
    ("minSNR", "minimum SNR value for loop output",
     cxxopts::value<double>(minSNR)->default_value("0.01"))
    ("prefix", "Filename prefix",
     cxxopts::value<string>(prefix)->default_value("crossval"))
    ("runtag", "Phase space file",
     cxxopts::value<string>(runtag)->default_value("run1"))
    ("outdir", "Output directory path",
     cxxopts::value<string>(outdir)->default_value("."))
    ("modelfile", "Halo model file",
     cxxopts::value<string>(modelf)->default_value("SLGridSph.model"))
    ("init", "fiducial phase-space index",
     cxxopts::value<int>(init)->default_value("0"))
    ("num", "Number of entries in Q table",
     cxxopts::value<int>(num)->default_value("10000"))
    ("compname", "train on Component (default=stars)",
     cxxopts::value<std::string>(cname)->default_value("stars"))
    ("psfile", "List of phase space files for processing",
     cxxopts::value<std::string>(psfiles))
    ("delimiter", "Phase-space file list delimiter for node index",
     cxxopts::value<std::string>(delim))
    ("ignore", "rebuild EOF grid if input parameters do not match the cachefile",
     cxxopts::value<bool>(ignore)->default_value("false"))
    ("cachefile", "cachefile name",
     cxxopts::value<std::string>(CACHEFILE)->default_value(".eof.cache.file"))
    ("tablefile", "table file name",
     cxxopts::value<std::string>(table_cache)->default_value(".cross_val_cyl"))
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
  if (vm.count("LOG")) LOG = true;

  bool Hall = false;
  if (vm.count("Hall")) Hall = true;

  bool verbose = false;
  if (vm.count("verbose")) verbose = true;

  // ==================================================
  // Nice process
  // ==================================================

  if (NICE>0)
    setpriority(PRIO_PROCESS, 0, NICE);

  std::ofstream out(prefix+".summary");
  if (not out) {
    std::cerr << "Error opening output file <" << prefix+".summary" << ">" << std::endl;
    MPI_Finalize();
    exit(-2);
  }

  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  int mmax, numx, numy, norder, cmapr, cmapz, nodd=-1;
  double rcylmin, rcylmax, vscale;

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
	if (node["odd"])
	  nodd  = node["odd"  ].as<int>();
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
      }
    }
  }

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAPR       = cmapr;
  EmpCylSL::CMAPZ       = cmapz;
  EmpCylSL::logarithmic = LOG;
  EmpCylSL::PCAVAR      = true;

				// Create expansion
				//
  EmpCylSL ortho(NMAX, LMAX, mmax, norder, rscale, vscale, nodd, CACHEFILE);
    
				// Set smoothing type to Truncate or
				// Hall (default)
  if (vm.count("truncate"))
    ortho.setTK("Truncate");
  else
    ortho.setTK("Hall");

  if (NPART) ortho.setSampT(NPART);

  std::vector<Particle> particles;
  PR::PRptr reader;
  
  std::vector<double> times;
  std::vector<std::string> outfiles;

  if (ortho.read_cache()==0) {
    std::cout << "Could not read cache file <" << CACHEFILE << ">"
	      << " . . . quitting" << std::endl;
    MPI_Finalize();
    exit(0);
  }
  

  // ==================================================
  // Compute and table Q values
  // ==================================================

  double Ascale = ortho.get_ascale();
  double Rtable = M_SQRT1_2 * RMAX;

  double XMIN = ortho.r_to_xi(RMIN*Ascale);
  double XMAX = ortho.r_to_xi(Rtable*Ascale);
  double dX   = (XMAX - XMIN)/EmpCylSL::NUMX;
    
  double YMIN = ortho.z_to_y(-Rtable*Ascale);
  double YMAX = ortho.z_to_y( Rtable*Ascale);
  double dY   = (YMAX - YMIN)/EmpCylSL::NUMY;
    
  // ============================
  // Begin array integration loop
  // ============================

  std::vector<Eigen::MatrixXd> Ec((mmax+1)*norder), Es((mmax+1)*norder);

  for (int M=0; M<=mmax; M++) {
    for (int n=0; n<norder; n++) {
      Ec[M*norder+n].resize(EmpCylSL::NUMX, EmpCylSL::NUMY+1);
      if (M) Es[M*norder+n].resize(EmpCylSL::NUMX, EmpCylSL::NUMY+1);
    }
  }


  std::ifstream in(table_cache);
  bool MakeCache = true;
    
  if (in) {

    const unsigned int cmagic = 0xf00ba;

    // Attempt to read magic number
    //
    unsigned int tmagic;
    in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));


    if (tmagic == cmagic) {

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
	if (myid==0)
	  std::cerr << "YAML: error parsing <" << buf.get() << "> "
		    << "in " << __FILE__ << ":" << __LINE__ << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
      }
      
      // Get parameters
      //
      int    MMAX   = node["mmax"  ].as<int>();
      int    NORDER = node["norder"].as<int>();
      int    numx   = node["numx"  ].as<int>();
      int    numy   = node["numy"  ].as<int>();
      double rmin   = node["rmin"  ].as<double>();
      double rmax   = node["rmax"  ].as<double>();
      double ascl   = node["ascl"  ].as<double>();

      bool okay = true;
      if (MMAX   != mmax          )   okay = false;
      if (NORDER != norder        )   okay = false;
      if (numx   != EmpCylSL::NUMX)   okay = false;
      if (numy   != EmpCylSL::NUMY)   okay = false;
      if (fabs(rmin-RMIN)   > 1.0e-8) okay = false;
      if (fabs(rmax-RMAX)   > 1.0e-8) okay = false;
      if (fabs(ascl-Ascale) > 1.0e-8) okay = false;

      // Read table
      //
      if (okay) {
	MakeCache = false;

	for (int M=0; M<=mmax; M++) {
	  for (int n=0; n<norder; n++) {
	    int id = M*norder + n;
	    in.read(reinterpret_cast<char *>(Es[id].data()), Es[id].size()*sizeof(double));
	    if (M)
	      in.read(reinterpret_cast<char *>(Es[id].data()), Es[id].size()*sizeof(double));
	  }
	}
      } 
    }
  }
  

  // Get data from EmpCylSL
  //
  auto densC = ortho.getDensC();
  auto densS = ortho.getDensS();
  auto potlC = ortho.getPotlC();
  auto potlS = ortho.getPotlS();

  if (MakeCache) {


    // Storage temp
    //
    std::vector<double> PL(mmax+2), QL(mmax+2);

    // Double sum on grid
    //
    int icnt = 0;
    for (int i=0; i<EmpCylSL::NUMX; i++) {
      for (int j=0; j<=EmpCylSL::NUMY; j++) {
	
	if (icnt++ % numprocs == myid) {
	  
	  double R = ortho.xi_to_r(XMIN + dX*(0.5+i));
	  double z = ortho. y_to_z(YMIN + dY*j);
	  
#ifdef DEBUG
	  std::cout << std::setw(4)  << myid
		    << std::setw(4)  << i
		    << std::setw(4)  << j
		    << std::setw(16) << R
		    << std::setw(16) << z
		    << std::endl;
#endif
	  
	  for (int k=0; k<=EmpCylSL::NUMX; k++) {
	    for (int l=0; l<=EmpCylSL::NUMY; l++) {
	      
	      double Rp = ortho.xi_to_r(XMIN + dX*k);
	      double zp = ortho. y_to_z(YMIN + dY*l);
	      double xi = Xi(R, Rp, z, zp);

	      if (not std::isinf(xi)) {

		int zero = 0, nord = mmax+1, newn;
		dtorh1_(&xi, &zero, &nord, &PL[0], &QL[0], &newn);
		
		for (int M=0; M<=mmax; M++) {
		  for (int n=0; n<norder; n++) {
		    int id = M*norder + n;
		    Ec[id](i, j) += -2.0*sqrt(Rp/R)*QL[M]*densC[id](k, l);
		    if (M)
		      Es[id](i, j) += -2.0*sqrt(Rp/R)*QL[M]*densS[id](k, l);
		  }
		}
	      }
	      // END: is not inf
	    }
	    // END: inner y loop
	  }
	  // END: inner x loop
	}
	// END: MPI node selection
      }
      // END: outer y loop
    }
    // END: outer x loop
  }
  // END: MakeCache

#ifdef DEBUG
  std::cout << "[" << myid << "] Before synchronization" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // MPI: share data
  //
  for (int M=0; M<=mmax; M++) {
    for (int n=0; n<norder; n++) {
      int id = M*norder + n;
      MPI_Allreduce(MPI_IN_PLACE, Ec[id].data(), Ec[id].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (M)
	MPI_Allreduce(MPI_IN_PLACE, Es[id].data(), Es[id].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }
  
#ifdef DEBUG
  std::cout << "[" << myid << "] After synchronization" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (myid==0 and MakeCache) {

    std::ofstream out(table_cache);
    
    if (out) {
      
#ifdef DEBUG
      std::cout << "[" << myid << "] Begin to write table cache" << std::endl;
#endif

      const unsigned int cmagic = 0xf00ba;

      // This is a node of simple {key: value} pairs.  More general
      // content can be added as needed.
      YAML::Node node;
    
      node["mmax"  ] = mmax;
      node["norder"] = norder;
      node["numx"  ] = EmpCylSL::NUMX;
      node["numy"  ] = EmpCylSL::NUMY;
      node["rmin"  ] = RMIN;
      node["rmax"  ] = RMAX;
      node["ascl"  ] = Ascale;
    
      // Serialize the node
      //
      YAML::Emitter y; y << node;
      
      // Get the size of the string
      //
      unsigned int hsize = strlen(y.c_str());
      
      // Write magic #
      //
      out.write(reinterpret_cast<const char *>(&cmagic),   sizeof(unsigned int));
      
      // Write YAML string size
      //
      out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));
      
      // Write YAML string
      //
      out.write(reinterpret_cast<const char *>(y.c_str()), hsize);
      
      
      // Write table
      //
      for (int M=0; M<=mmax; M++) {
	for (int n=0; n<norder; n++) {
	  int id = M*norder + n;
	  out.write(reinterpret_cast<const char *>(Es[id].data()), Es[id].size()*sizeof(double));
	  if (M)
	    out.write(reinterpret_cast<const char *>(Es[id].data()), Es[id].size()*sizeof(double));
	}
      }

#ifdef DEBUG
      std::cout << "[" << myid << "] Wrote cache file" << std::endl;
#endif

    } else {
      std::cout << "Could not open table cache <" << table_cache << ">" << std::endl;
    }

#ifdef DEBUG
    std::cout << "[" << myid << "] Finished writing table cache" << std::endl;
#endif
  }

  // ==================================================
  // Phase space output loop
  // ==================================================

  std::string file;

#ifdef DEBUG
  std::cout << "[" << myid << "] Begin phase-space loop" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif	      

  unsigned ibatch = 0;

  for (auto batch : PR::ParticleReader::parseFileList(psfiles, delim)) {
    
    PR::PRptr reader = PR::ParticleReader::createReader
      (fileType, batch, myid, true);

    double tnow = reader->CurrentTime();
    if (myid==0) std::cout << "Beginning partition [time=" << tnow
			   << ", index=" << ibatch << "] . . . "  << flush;
    
    reader->SelectType(cname);
      
    //------------------------------------------------------------ 

    if (myid==0) std::cout << std::endl
			   << "Accumulating particle positions . . . "
			   << std::flush;

    ortho.setup_accumulation();
    ortho.setHall("test", reader->CurrentNumber());

    auto p = reader->firstParticle();
    int icnt = 0;
    do {
      if (icnt++ % numprocs == myid) {
	double R   = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	double phi = atan2(p->pos[1], p->pos[0]);
	ortho.accumulate(R, p->pos[2], phi, p->mass, p->indx, 0, 0, true);
	//                                                    ^  ^  ^
	//                                                    |  |  |
	// Thread id -----------------------------------------+  |  |
	// Level ------------------------------------------------+  |
	// Compute covariance --------------------------------------+
      }
      p = reader->nextParticle();
    } while (p);
    
    if (myid==0) std::cout << "done" << endl;
    
    //------------------------------------------------------------ 
      
    if (myid==0) cout << "Making coefficients . . . " << flush;
    ortho.make_coefficients(true);
    ortho.pca_hall(true, false);
    if (myid==0) std::cout << "done" << endl;

    //------------------------------------------------------------ 
    

    std::vector<double> term1(mmax+1);
    std::vector<double> term2(mmax+1), work2(mmax+1);
    std::vector<double> term3(mmax+1), work3(mmax+1);
    
    double minSNR = ortho.getMinSNR();
    double maxSNR = ortho.getMaxSNR();

    if (myid==0) {
      std::cout << "Found minSNR=" << minSNR
		<< " maxSNR=" << maxSNR << std::endl;
    }

    if (maxSNR < minSNR)        minSNR = maxSNR * 1.0e-2;
    if (minSNR < maxSNR*1.0e-6) minSNR = maxSNR * 1.0e-6;
    
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

    double term4tot = 0.0;
    std::vector<Eigen::VectorXd> ac_cos, ac_sin;

    for (int nsnr=0; nsnr<NSNR; nsnr++) {

      // Assign the snr value
      //
      double snr = minSNR + dSNR*nsnr;
      if (LOG) snr = exp(snr);

      if (myid==0) {
	std::cout << "Computing SNR=" << snr;
	if (Hall) std::cout << " using Hall smoothing . . . " << std::endl;
	else      std::cout << " using truncation . . . " << std::endl;
      }
    
      // Get the snr trimmed coefficients
      //
      ortho.get_trimmed(snr, ac_cos, ac_sin);
	
      if (myid==0) {
	std::cout << "*** SNR: " << snr << std::endl;
	for (int M=0; M<=mmax; M++) {
	  for (int i=0; i<norder; i++)
	    std::cout << std::setw(4)  << M
		      << std::setw(4)  << i
		      << std::setw(18) << ac_cos[M][i] << std::endl;
	  std::cout << std::endl;
	}
      }

      // Zero out the accumulators
      //
      std::fill(term1.begin(), term1.end(), 0.0);
      std::fill(term2.begin(), term2.end(), 0.0);
      std::fill(term3.begin(), term3.end(), 0.0);
      std::fill(work2.begin(), work2.end(), 0.0);
      std::fill(work3.begin(), work3.end(), 0.0);

      if (myid==0) {		// Only root process needs this one
	  
	// Term 1
	//
	double t1 = 0.0;
	for (int M=0; M<=mmax; M++) {
	  for (int n=0; n<norder; n++) {
	    if (M==0)
	      term1[M] += ac_cos[M][n]*ac_cos[M][n];
	    else {
	      term1[M] +=
		ac_cos[M][n]*ac_cos[M][n] + 
		ac_sin[M][n]*ac_sin[M][n] ;
	    }
	  }
	}
      }

      // Term 2 and Term 3
      //
      for (int M=0; M<=mmax; M++) {

	// Particle loop
	//
	p = reader->firstParticle();
	int icnt = 0;
	do {
	  if (icnt++ % numprocs == myid) {
	    double r = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	    double z = p->pos[2];
	    double phi = atan2(p->pos[1], p->pos[0]);
	    double mass = p->mass;
	    
	    double x = ortho.r_to_xi(r);
	    x = std::max<double>(XMIN, x);
	    x = std::min<double>(XMAX, x);
	    
	    int iX = floor( (x - XMIN - 0.5*dX)/dX );
	    if (iX<0) iX = 0;
	    if (iX>=EmpCylSL::NUMX-1) iX = EmpCylSL::NUMX-2;
	    
	    double y = ortho.z_to_y(z);
	    y = std::max<double>(YMIN, y);
	    y = std::min<double>(YMAX, y);
	    
	    int iY = floor( (y - YMIN)/dY );
	    if (iY<0) iY = 0;
	    if (iY>=EmpCylSL::NUMY) iY = EmpCylSL::NUMY-1;
	    
#ifdef DEBUG
           if (iX<0 or iX>=EmpCylSL::NUMX-1) {
              std::cout << "X out of bounds: x=" << x << " iX=" << iX
                        << " XMIN=" << XMIN
                        << " XMAX=" << XMAX
                        << std::endl;
            }
	   
            if (iY<0 or iY>=EmpCylSL::NUMY) {
              std::cout << "Y out of bounds: y=" << y << " iY=" << iY
                        << " YMIN=" << YMIN
                        << " YMAX=" << YMAX
                        << std::endl;
            }
#endif
				// X dimension
	    double A = (XMIN + dX*(iX+1.5) - x)/dX;
	    double B = (x - XMIN - dX*(iX+0.5))/dX;
				// Y dimension
	    double C = (YMIN + dY*(iY+1) - y)/dY;
	    double D = (y - YMIN - dY*(iY+0))/dY;
	    
	    for (int M=0; M<=mmax; M++) {

	      for (int n=0; n<norder; n++) {

		int id = M*norder + n;
		double PotlS = 0.0, DensS = 0.0;

		double PotlC =
		  A*(C*potlC[id](iX+0, iY+0) + D*potlC[id](iX+0, iY+1) ) +
		  B*(C*potlC[id](iX+1, iY+0) + D*potlC[id](iX+1, iY+1) ) ;

		double DensC =
		  A*(C*Ec[id](iX+0, iY+0) + D*Ec[id](iX+0, iY+1) ) +
		  B*(C*Ec[id](iX+1, iY+0) + D*Ec[id](iX+1, iY+1) ) ;

		if (M) {
		  PotlS =
		    A*(C*potlS[id](iX+0, iY+0) + D*potlS[id](iX+0, iY+1) ) +
		    B*(C*potlS[id](iX+1, iY+0) + D*potlS[id](iX+1, iY+1) ) ;

		  DensS =
		    A*(C*Es[id](iX+0, iY+0) + D*Es[id](iX+0, iY+1) ) +
		    B*(C*Es[id](iX+1, iY+0) + D*Es[id](iX+1, iY+1) ) ;
		}

		if (M==0) {
		  work2[M] += mass*ac_cos[M][n]*PotlC;
		  work3[M] += mass*ac_cos[M][n]*DensC;
		} else {
		  work2[M] += mass*(ac_cos[M][n]*PotlC*cos(phi*M) + ac_sin[M][n]*PotlS*sin(phi*M));
		  work3[M] += mass*(ac_cos[M][n]*DensC*cos(phi*M) + ac_sin[M][n]*DensS*sin(phi*M));
		}
	      }
	      // END: radial index loop
	    }
	    // END: M loop
	  }
	  // END: add particle data
	    
	  // Queue up next particle
	  //
	  p = reader->nextParticle();
	} while (p);
	//
	// END: particle loop
      }
      // END: L loop

	
      MPI_Reduce(work2.data(), term2.data(), work2.size(), MPI_DOUBLE,
		 MPI_SUM, 0, MPI_COMM_WORLD);

      MPI_Reduce(work3.data(), term3.data(), work3.size(), MPI_DOUBLE,
		 MPI_SUM, 0, MPI_COMM_WORLD);

      if (myid==0) {
	  
	constexpr double pi4 = 4.0*M_PI;

	out << std::setw( 5) << ibatch
	    << std::setw(18) << snr;
	
	double term1tot = std::accumulate(term1.begin(), term1.end(), 0.0) / pi4;
	double term2tot = std::accumulate(term2.begin(), term2.end(), 0.0) * (-1);
	double term3tot = std::accumulate(term3.begin(), term3.end(), 0.0) ;

	if (nsnr==0) term4tot = term1tot;
	  
	out << std::setw(18) << term1tot
	    << std::setw(18) << term2tot
	    << std::setw(18) << term3tot
	    << std::setw(18) << term1tot - term2tot - term3tot + term4tot
	    << std::endl;
      }
      // Root process
      
      if (myid==0) std::cout << "done" << endl;
    }
    // SNR loop
      
    // Blank line between stanzas
    //
    if (myid==0) out << std::endl;

    ibatch++;
    
  } // Dump loop

  MPI_Finalize();

  return 0;
}

