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
#include <map>

using namespace std;

				// Eigen3
#include <Eigen/Eigen>

				// Boost stuff

#include <boost/shared_ptr.hpp>
#include <boost/make_unique.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

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
  
// Globals
//

extern double Ylm01(int ll, int mm);
extern double plgndr(int, int, double);

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
  int NICE, LMAX, NMAX, NSNR, NPART;
  int beg, end, stride, init, knots, num;
  std::string CACHEFILE, modelf, dir("./"), cname, prefix, table_cache;
  bool ignore;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Cross-validation analysis for cylindrical models" << std::endl
       << std::string(60, '-') << std::endl << std::endl
       << "Allowed options";
  
  po::options_description desc(sout.str());
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("verbose,v",                                                                       "Verbose and diagnostic output for covariance computation")
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
    ("RMIN",                po::value<double>(&RMIN)->default_value(0.0),
     "minimum radius for output")
    ("RSCALE",              po::value<double>(&rscale)->default_value(0.067),
     "coordinate mapping scale factor")
    ("RMAX",                po::value<double>(&RMAX)->default_value(2.0),
     "maximum radius for output")
    ("LMAX",                po::value<int>(&LMAX)->default_value(4),
     "Maximum harmonic order for spherical expansion")
    ("NMAX",                po::value<int>(&NMAX)->default_value(12),
     "Maximum radial order for spherical expansion")
    ("NPART",               po::value<int>(&NPART)->default_value(0),
     "Jackknife partition number for testing (0 means off, use standard eval)")
    ("NSNR, N",             po::value<int>(&NSNR)->default_value(20),
     "Number of SNR evaluations")
    ("minSNR",              po::value<double>(&minSNR)->default_value(0.01),
     "minimum SNR value for loop output")
    ("prefix",              po::value<string>(&prefix)->default_value("crossval"),
     "Filename prefix")
    ("runtag",              po::value<string>(&runtag)->default_value("run1"),
     "Phase space file")
    ("outdir",              po::value<string>(&outdir)->default_value("."),
     "Output directory path")
    ("modelfile",           po::value<string>(&modelf)->default_value("SLGridSph.model"),
     "Halo model file")
    ("init",                po::value<int>(&init)->default_value(0),
     "fiducial PSP index")
    ("beg",                 po::value<int>(&beg)->default_value(0),
     "initial PSP index")
    ("end",                 po::value<int>(&end)->default_value(99999),
     "final PSP index")
    ("stride",              po::value<int>(&stride)->default_value(1),
     "PSP index stride")
    ("num",                 po::value<int>(&num)->default_value(10000),
     "Number of entries in Q table")
    ("knots",               po::value<int>(&knots)->default_value(40),
     "Number of Legendre integration knots")
    ("compname",            po::value<std::string>(&cname)->default_value("stars"),
     "train on Component (default=stars)")
    ("dir,d",               po::value<std::string>(&dir),
     "directory for SPL files")
    ("ignore",
     po::value<bool>(&ignore)->default_value(false),
     "rebuild EOF grid if input parameters do not match the cachefile")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("tablefile",
     po::value<std::string>(&table_cache)->default_value(".cross_val_cyl"),
     "table file name")
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
    exit(-1);
  }

  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    if (myid==0) std::cout << std::endl << desc << std::endl;
    return 0;
  }

  bool SPL = false;
  if (vm.count("SPL")) SPL = true;
  if (vm.count("OUT")) SPL = false;

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
    exit(-2);
  }

  // ==================================================
  // PSP input stream
  // ==================================================

  int iok = 1;
  ifstream in0, in1;
  std::ostringstream s0, s1;
  if (myid==0) {
    s0 << "OUT." << runtag << "."
       << std::setw(5) << std::setfill('0') << init;
    in0.open(s0.str());
    if (!in0) {
      cerr << "Error opening <" << s0.str() << ">" << endl;
      iok = 0;
    }
  }
    
  MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (iok==0) {
    MPI_Finalize();
    exit(-1);
  }
    
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
	in.read((char *)&NMAX,    sizeof(int));
	in.read((char *)&norder,  sizeof(int));
	in.read((char *)&DENS,    sizeof(int)); 
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
  EmpCylSL::logarithmic = logl;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::CACHEFILE   = CACHEFILE;
  EmpCylSL::PCAVAR      = true;

				// Create expansion
				//
  EmpCylSL ortho(NMAX, LMAX, mmax, norder, rscale, vscale);
    
				// Set smoothing type to truncate
				//
  ortho.setTK("Truncate");


  vector<Particle> particles;
  PSPptr psp;
  
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

  std::map< std::pair<int, int>, Eigen::MatrixXd > Ec, Es;

  for (int M=0; M<=mmax; M++) {
    for (int n=0; n<norder; n++) {
      std::pair<int, int> id(M, n);
      Ec[id].resize(EmpCylSL::NUMX, EmpCylSL::NUMY+1);
      if (M) Es[id].resize(EmpCylSL::NUMX, EmpCylSL::NUMY+1);
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
      auto buf = boost::make_unique<char[]>(ssize+1);
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
	    std::pair<int, int> id(M, n);
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
    std::vector<double> PL(mmax+1), QL(mmax+1);

    // Double sum on grid
    //
    int icnt = 0;
    for (int i=0; i<EmpCylSL::NUMX; i++) {
      for (int j=0; j<=EmpCylSL::NUMY; j++) {
	
	if (icnt++ % numprocs == myid) {
	  
	  double R = ortho.xi_to_r(XMIN + dX*(0.5+i));
	  double z = ortho. y_to_z(YMIN + dY*j);
	  
	  std::cout << std::setw(4)  << myid
		    << std::setw(4)  << i
		    << std::setw(4)  << j
		    << std::setw(16) << R
		    << std::setw(16) << z
		    << std::endl;
	  
	  for (int k=0; k<=EmpCylSL::NUMX; k++) {
	    for (int l=0; l<=EmpCylSL::NUMY; l++) {
	      
	      double Rp = ortho.xi_to_r(XMIN + dX*k);
	      double zp = ortho. y_to_z(YMIN + dY*l);
	      
	      /*
	      std::cout << std::setw(44) << "***"
			<< std::setw(4)  << myid
			<< std::setw(4)  << l
			<< std::setw(4)  << k
			<< std::setw(16) << Rp
			<< std::setw(16) << zp
			<< std::endl;
	      */

	      double xi = Xi(R, Rp, z, zp);

	      if (not std::isinf(xi)) {

		int zero = 0, nord = mmax+1, newn;
		dtorh1_(&xi, &zero, &nord, &PL[0], &QL[0], &newn);
		
		for (int M=0; M<=mmax; M++) {
		  for (int n=0; n<norder; n++) {
		    std::pair<int, int> id(M, n);
		    Ec[id](i, j) += -2.0*sqrt(Rp/R)*QL[M]*densC[id](k, l);
		    if (M)
		      Es[id](i, j) += -2.0*sqrt(Rp/R)*QL[M]*densS[id](k, l);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  std::cout << "[" << myid << "] Before synchronization" << std::endl;

  // MPI share data
  //
  for (int M=0; M<=mmax; M++) {
    for (int n=0; n<norder; n++) {
      std::pair<int, int> id(M, n);
      MPI_Allreduce(MPI_IN_PLACE, Ec[id].data(), Ec[id].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (M)
	MPI_Allreduce(MPI_IN_PLACE, Es[id].data(), Es[id].size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }
  
  std::cout << "[" << myid << "] After synchronization" << std::endl;

  if (myid==0 and MakeCache) {

    std::ofstream out(table_cache);
    
    if (out) {

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
      out.write(reinterpret_cast<const char *>(cmagic),   sizeof(unsigned int));
      
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
	  std::pair<int, int> id(M, n);
	  out.write(reinterpret_cast<const char *>(Es[id].data()), Es[id].size()*sizeof(double));
	  if (M)
	    out.write(reinterpret_cast<const char *>(Es[id].data()), Es[id].size()*sizeof(double));
	}
      }

      std::cout << "[" << myid << "] Wrote cache file" << std::endl;

    } else {
      std::cout << "Could not open table cache <" << table_cache << ">" << std::endl;
    }
  }

  // ==================================================
  // Phase space output loop
  // ==================================================

  std::string file;

  std::cout << "[" << myid << "] Begin phase-space loop" << std::endl;


  for (int ipsp=beg; ipsp<=end; ipsp+=stride) {

    // ==================================================
    // PSP input stream
    // ==================================================

    int iok = 1;
    std::ostringstream s0, s1;

    s1.str("");		// Clear stringstream
    if (SPL) s1 << "SPL.";
    else     s1 << "OUT.";
    s1 << runtag << "."<< std::setw(5) << std::setfill('0') << ipsp;
      
				// Check for existence of next file
    file = dir + s1.str();
    std::ifstream in(file);
    if (!in) {
      std::cerr << "Error opening <" << file << ">" << endl;
      iok = 0;
    }
    
    if (iok==0) break;

    // ==================================================
    // Open PSP file
    // ==================================================
    PSPptr psp;

    if (SPL) psp = std::make_shared<PSPspl>(s1.str(), dir, true);
    else     psp = std::make_shared<PSPout>(file, true);

    tnow = psp->CurrentTime();
    if (myid==0) std::cout << "Beginning partition [time=" << tnow
			   << ", index=" << ipsp << "] . . . "  << flush;
    
    if (not psp->GetNamed(cname)) {
      if (myid==0) {
	std::cout << "Error finding component named <" << cname << ">" << std::endl;
	psp->PrintSummary(std::cout);
      }
      exit(-1);
    }
      
    //------------------------------------------------------------ 

    if (myid==0) std::cout << std::endl
			   << "Accumulating particle positions . . . "
			   << std::flush;

    ortho.setup_accumulation();

    SParticle *p = psp->GetParticle();
    int icnt = 0;
    do {
      if (icnt++ % numprocs == myid) {
	double R   = sqrt(p->pos(0)*p->pos(0) + p->pos(0)*p->pos(1));
	double phi = atan2(p->pos(1), p->pos(0));
	ortho.accumulate(R, p->pos(2), phi, p->mass(), p->indx(), 0);
      }
      p = psp->NextParticle();
    } while (p);
    
    if (myid==0) std::cout << "done" << endl;
    
    //------------------------------------------------------------ 
      
    if (myid==0) cout << "Making coefficients . . . " << flush;
    ortho.make_coefficients();
    if (myid==0) std::cout << "done" << endl;

    //------------------------------------------------------------ 
    

    std::vector<double> term1(mmax+1);
    std::vector<double> term2(mmax+1), work2(mmax+1);
    std::vector<double> term3(mmax+1), work3(mmax+1);
    
    double maxSNR = ortho.getMaxSNR();

    if (maxSNR < minSNR) minSNR = maxSNR / 100.0;
    
    if (LOG) {
      minSNR = log(minSNR);
      maxSNR = log(maxSNR);
    }

    double dSNR = (maxSNR - minSNR)/(NSNR - 1);

    if (myid==0) {
      std::cout << "maxSNR=" << maxSNR << " dSNR=" << dSNR << std::endl;
    }

    double term4tot = 0.0;

    for (int nsnr=0; nsnr<NSNR; nsnr++) {

      // Assign the snr value
      //
      double snr = minSNR + dSNR*nsnr;
      if (LOG) snr = exp(snr);

      if (myid==0) {
	std::cout << "Computing SNR=" << snr;
	if (Hall) std::cout << " using Hall smoothing . . . " << flush;
	else      std::cout << " using truncation . . . " << flush;
      }
    
      // Get the snr trimmed coefficients
      //
      std::vector<Vector> ac_cos, ac_sin;
      ortho.get_trimmed(snr, ac_cos, ac_sin);
	
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
	p = psp->GetParticle();
	int icnt = 0;
	do {
	  if (icnt++ % numprocs == myid) {
	    double r = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
	    double z = p->pos(2);
	    double phi = atan2(p->pos(1), p->pos(0));
	    double mass = p->mass();
	    
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
	    
	    double A = (XMIN + dX*(iX+1.5) - x)/dX;
	    double B = (x - XMIN - dX*(iX+0.5))/dX;
	    
	    double C = (YMIN + dY*(iY+1) - y)/dY;
	    double D = (y - YMIN - dY*(iY+0))/dY;
	    
	    for (int M=0; M<=mmax; M++) {

	      for (int n=0; n<norder; n++) {

		std::pair<int, int> id(M, n);
		double PotlS = 0.0, DensS = 0.0;

		double PotlC =
		  A*(C*potlC[id](iX+0, iY+0) + D*potlC[id](iX+0, iY+1) ) +
		  B*(C*potlC[id](iX+1, iY+0) + D*potlC[id](iX+1, iY+1) ) ;

		double DensC =
		  A*(C*densC[id](iX+0, iY+0) + D*densC[id](iX+0, iY+1) ) +
		  B*(C*densC[id](iX+1, iY+0) + D*densC[id](iX+1, iY+1) ) ;

		if (M) {
		  PotlS =
		    A*(C*potlS[id](iX+0, iY+0) + D*potlS[id](iX+0, iY+1) ) +
		    B*(C*potlS[id](iX+1, iY+0) + D*potlS[id](iX+1, iY+1) ) ;

		  DensS =
		    A*(C*densS[id](iX+0, iY+0) + D*densS[id](iX+0, iY+1) ) +
		    B*(C*densS[id](iX+1, iY+0) + D*densS[id](iX+1, iY+1) ) ;
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
	  p = psp->NextParticle();
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
	  
	out << std::setw( 5) << ipsp
	    << std::setw(18) << snr;
	
	double term1tot = std::accumulate(term1.begin(), term1.end(), 0.0) / (4.0*M_PI);
	double term2tot = std::accumulate(term2.begin(), term2.end(), 0.0);
	double term3tot = std::accumulate(term3.begin(), term3.end(), 0.0);

	if (nsnr==0) term4tot = term1tot;
	  
	out << std::setw(18) << term1tot
	    << std::setw(18) << term2tot
	    << std::setw(18) << term3tot
	    << std::setw(18) << term1tot + term2tot - term3tot + term4tot
	    << std::endl;
      }
      // Root process
      
      if (myid==0) std::cout << "done" << endl;

    }
    // SNR loop
      
    // Blank line between stanzas
    //
    if (myid==0) out << std::endl;

  } // Dump loop

  MPI_Finalize();

  return 0;
}

