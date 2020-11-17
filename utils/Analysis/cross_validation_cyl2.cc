/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Cross validation analysis for cylinder
 *
 *  Monopole expansion version
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
  int beg, end, stride, init, knots, numr;
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
    ("LMAX",                po::value<int>(&LMAX)->default_value(24),
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
    ("numr",                po::value<int>(&numr)->default_value(256),
     "Number of entries in R table")
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
  double dX0  = (XMAX - XMIN)/EmpCylSL::NUMX;
    
  double YMIN = ortho.z_to_y(-Rtable*Ascale);
  double YMAX = ortho.z_to_y( Rtable*Ascale);
  double dY0  = (YMAX - YMIN)/EmpCylSL::NUMY;

  double dX   = (XMAX - XMIN)/(numr+1);
    
  // ============================
  // Begin array integration loop
  // ============================

  std::vector<std::vector<double>> Ec((LMAX+1)*(mmax+1)*norder);
  std::vector<std::vector<double>> Es((LMAX+1)*(mmax+1)*norder);

  for (int L=0; L<=LMAX; L++) {
    for (int M=0; M<=std::min<int>(L, mmax); M++) {
      for (int n=0; n<norder; n++) {
	int id = (L*(mmax+1) + M)*norder + n;
	Ec[id].resize(numr+1);
	if (M) Es[id].resize(numr+1);
      }
    }
  }

  std::ifstream in(table_cache);
  bool MakeCache = true;
    
  if (in) {

    const unsigned int cmagic = 0xf00bb;

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
      int    lmax   = node["lmax"  ].as<int>();
      int    MMAX   = node["mmax"  ].as<int>();
      int    NORDER = node["norder"].as<int>();
      int    NUMR   = node["numr"  ].as<int>();
      double rmin   = node["rmin"  ].as<double>();
      double rmax   = node["rmax"  ].as<double>();
      double ascl   = node["ascl"  ].as<double>();

      bool okay = true;
      if (LMAX   != lmax          )   okay = false;
      if (MMAX   != mmax          )   okay = false;
      if (NORDER != norder        )   okay = false;
      if (NUMR   != numr          )   okay = false;
      if (fabs(rmin-RMIN)   > 1.0e-8) okay = false;
      if (fabs(rmax-RMAX)   > 1.0e-8) okay = false;
      if (fabs(ascl-Ascale) > 1.0e-8) okay = false;

      // Read table
      //
      if (okay) {
	MakeCache = false;

	for (int L=0; L<=LMAX; L++) {
	  for (int M=0; M<=std::min<int>(L, mmax); M++) {
	    for (int n=0; n<norder; n++) {
	      int id = (L*(mmax+1) + M)*norder + n;
	      in.read(reinterpret_cast<char *>(Es[id].data()), Es[id].size()*sizeof(double));
	      if (M) 
		in.read(reinterpret_cast<char *>(Es[id].data()), Es[id].size()*sizeof(double));
	    }
	  }
	}
      } 
    }
  }
  

  // Get data from EmpCylSL
  //
  auto potlC = ortho.getPotlC();
  auto potlS = ortho.getPotlS();

  if (MakeCache) {

    const int NINTR = 40;
    const int NINTT = 40;

    LegeQuad lr(NINTR);
    LegeQuad lt(NINTT);

    // Inner
    //
    int icnt = 0;

    for (int L=0; L<=LMAX; L++) {

      for (int M=0; M<=std::min<int>(L, mmax); M++) {

	for (int n=0; n<norder; n++) {

	  int id = (L*(mmax+1) + M)*norder + n;

	  for (int i=0; i<=numr; i++) {
      
	    if (icnt++ % numprocs == myid) {
	  
	      double xi = XMIN + dX*i;
	      double ri = ortho.xi_to_r(xi);

	      double innerC = 0.0, outerC = 0.0;
	      double innerS = 0.0, outerS = 0.0;

	      for (int t=1; t<=NINTT; t++) {

		double cosx = -1.0 + 2.0*lt.knot(t);
		double sinx = sqrt(1.0 - cosx*cosx);
		double ylm  = Ylm01(L, M) * plgndr(L, M, cosx) *
		  4.0*M_PI/(2.0*L+1.0) * lt.weight(t);

		for (int k=1; k<=NINTR; k++) {
		  // Inner
		  //
		  double x = XMIN + (xi - XMIN)*lr.knot(k);
		  double r = ortho.xi_to_r(x);
		  double R = r*sinx, z = r*cosx;

		  double dC, dS;
		  ortho.getDensSC(M, n, r*sinx, r*cosx, dC, dS);

		  double fac = pow(r/ri, 1.0+L)*r * ylm *
		    (xi - XMIN)*lr.weight(k) / ortho.d_xi_to_r(x);
		
		  innerC += dC * fac;
		  innerS += dS * fac;

		  // Outer
		  //
		  x = xi + (XMAX - xi)*lr.knot(k);
		  r = ortho.xi_to_r(x);
		  R = r*sinx;
		  z = r*cosx;

		  ortho.getDensSC(M, n, r*sinx, r*cosx, dC, dS);

		  fac = pow(ri/r, L)*r * ylm *
		    (xi - XMIN)*lr.weight(k) / ortho.d_xi_to_r(x);
		
		  outerC += dC * fac;
		  outerS += dS * fac;
		}
	      }

	      Ec[id][i] += innerC + outerC;
	      if (M) Es[id][i] += innerS + outerS;
	    }
	    // END: MPI node selection
	  }
	  // END: radial grid
	}
	// END: n-order
      }
      // END: M value
    }
    // END: L value
  }
  // END: MakeCache

#ifdef DEBUG
  std::cout << "[" << myid << "] Before synchronization" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // MPI: share data
  //
  for (int L=0; L<=LMAX; L++) {
    for (int M=0; M<=std::min<int>(L, mmax); M++) {
      for (int n=0; n<norder; n++) {
	int id = (L*(mmax+1) + M)*norder + n;
	MPI_Allreduce(MPI_IN_PLACE, Ec[id].data(), Ec[id].size(), MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	if (M)
	  MPI_Allreduce(MPI_IN_PLACE, Es[id].data(), Es[id].size(), MPI_DOUBLE,
			MPI_SUM, MPI_COMM_WORLD);
      }
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

      const unsigned int cmagic = 0xf00bb;

      // This is a node of simple {key: value} pairs.  More general
      // content can be added as needed.
      YAML::Node node;
    
      node["lmax"  ] = LMAX;
      node["mmax"  ] = mmax;
      node["norder"] = norder;
      node["numr"  ] = numr;
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
      for (int L=0; L<=LMAX; L++) {
	for (int M=0; M<=std::min<int>(L, mmax); M++) {
	  for (int n=0; n<norder; n++) {
	    int id = (L*(mmax+1) + M)*norder + n;
	    out.write(reinterpret_cast<const char *>(Es[id].data()), Es[id].size()*sizeof(double));
	      if (M) 
		out.write(reinterpret_cast<const char *>(Es[id].data()), Es[id].size()*sizeof(double));
	  }
	}
      }
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
    ortho.setTotal(psp->GetNamed(cname)->comp.nbod); // Set particle number
    ortho.init_pca();

    SParticle *p = psp->GetParticle();
    int icnt = 0;
    do {
      if (icnt++ % numprocs == myid) {
	double R   = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
	double phi = atan2(p->pos(1), p->pos(0));
	ortho.accumulate(R, p->pos(2), phi, p->mass(), p->indx(), 0, 0, true);
	//                                                        ^  ^  ^
	//                                                        |  |  |
	// Thread id ---------------------------------------------+  |  |
	// Level ----------------------------------------------------+  |
	// Compute covariance ------------------------------------------+
      }
      p = psp->NextParticle();
    } while (p);
    
    if (myid==0) std::cout << "done" << endl;
    
    //------------------------------------------------------------ 
      
    if (myid==0) cout << "Making coefficients . . . " << flush;
    ortho.make_coefficients(true);
    ortho.pca_hall(true);
    if (myid==0) std::cout << "done" << endl;

    //------------------------------------------------------------ 
    

    std::vector<double> term1(mmax+1);
    std::vector<double> term2(mmax+1), work2(mmax+1);
    std::vector<double> term3(mmax+1), work3(mmax+1);
    
    double maxSNR = ortho.getMaxSNR();

    maxSNR = 5.0;

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
    std::vector<Vector> ac_cos, ac_sin;

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

      // Particle loop
      //
      p = psp->GetParticle();
      int icnt = 0;
      do {
	if (icnt++ % numprocs == myid) {
	  double R = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
	  double z = p->pos(2);
	  double phi = atan2(p->pos(1), p->pos(0));
	  double mass = p->mass();
	    
	  // Term 2
	  //
	  double x = ortho.r_to_xi(R);
	  x = std::max<double>(XMIN, x);
	  x = std::min<double>(XMAX, x);
	    
	  int iX = floor( (x - XMIN)/dX0 );
	  if (iX<0) iX = 0;
	  if (iX>=EmpCylSL::NUMX) iX = EmpCylSL::NUMX;
	    
	  double y = ortho.z_to_y(z);
	  y = std::max<double>(YMIN, y);
	  y = std::min<double>(YMAX, y);
	  
	  int iY = floor( (y - YMIN)/dY0 );
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
	  double A = (XMIN + dX0*(iX+1) - x)/dX0;
	  double B = (x - XMIN - dX0*(iX+0))/dX0;

	  // Y dimension
	  double C = (YMIN + dY0*(iY+1) - y)/dY0;
	  double D = (y - YMIN - dY0*(iY+0))/dY0;
	    
	  for (int M=0; M<=mmax; M++) {

	    for (int n=0; n<norder; n++) {

	      int id = M*norder + n;
	      double PotlS = 0.0;

	      double PotlC =
		A*(C*potlC[id](iX+0, iY+0) + D*potlC[id](iX+0, iY+1) ) +
		B*(C*potlC[id](iX+1, iY+0) + D*potlC[id](iX+1, iY+1) ) ;
	      
	      
	      if (M)
		PotlS =
		  A*(C*potlS[id](iX+0, iY+0) + D*potlS[id](iX+0, iY+1) ) +
		  B*(C*potlS[id](iX+1, iY+0) + D*potlS[id](iX+1, iY+1) ) ;
		
			if (M==0) 
		  work2[M] += mass*ac_cos[M][n]*PotlC;
		else 
		  work2[M] += mass*(ac_cos[M][n]*PotlC*cos(phi*M) + ac_sin[M][n]*PotlS*sin(phi*M));
		
	    }
	    // END: n order loop
	  }
	  // END: M loop

	  // Term 3
	  //
	  double r = sqrt(R*R + z*z);
	  double cosx = z/r;
	  x = ortho.r_to_xi(r);
	  x = std::max<double>(XMIN, x);
	  x = std::min<double>(XMAX, x);
	  
	  iX = floor( (x - XMIN)/dX );
	  if (iX<0) iX = 0;
	  if (iX>=numr) iX = numr - 1;
	  
	  A = (XMIN + dX*(iX+1) - x)/dX;
	  B = (x - XMIN - dX*(iX+0))/dX;
	    
	  for (int L=0; L<=LMAX; L++) {

	    for (int M=0; M<=std::min<int>(L, mmax); M++) {

	      double Ylm  = Ylm01(L, M) * plgndr(L, M, cosx);
	      double cosp = cos(phi*M), sinp = sin(phi*M);

	      for (int n=0; n<norder; n++) {

		int id = (L*(mmax+1) + M)*norder + n;
		
		work3[M] += mass * Ylm * ac_cos[M][n] * cosp * 
		  (A*Ec[id][iX+0] + B*Ec[id][iX+1]);
		if (M)
		  work3[M] += mass * Ylm * ac_sin[M][n] * sinp * 
		    (A*Es[id][iX+0] + B*Es[id][iX+1]);
	      }
	      // END: n order loop
	    }
	    // END: M loop
	  }
	  // END: L loop
	}
	// END: add particle data
	    
	// Queue up next particle
	//
	p = psp->NextParticle();
      } while (p);
      //
      // END: particle loop

	
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

