/*
  Separate a psp structure and make a kinematic Fourier coefficients
  series in Laguerre functions

  MDWeinberg 08/13/20
*/

#include <cstdlib>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>
#include <list>
#include <map>

#include <header.H>
#include <PSP.H>
#include <FileUtils.H>

#include <Progress.H>

#include <boost/program_options.hpp>
#include <boost/math/special_functions/laguerre.hpp>

#include <mpi.h>

namespace po = boost::program_options;

//
// MPI variables
//
int numprocs, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];



//! Generate orthonormal Laguerre functions
class Laguerre
{
private:

  double rscl;
  unsigned int nmax;
  std::vector<double> norm;

public:

  //! Constructor: set the order and number of radial functions
  Laguerre(double rscl, unsigned int nmax) : rscl(rscl), nmax(nmax)
  {
    norm.resize(nmax);
    for (unsigned int n=0; n<nmax; n++) norm[n] = 0.5*rscl*sqrt(1.0+n);
  }
  
  //! Get the norm for radial order m
  double getNorm(int n)
  {
    if (n>=nmax) return 0.0;
    else         return norm[n];
  }
  
  //! Evaluate the orthogonal Laguerre polynomial
  double operator()(double& r, const unsigned& n)
  {
    if (n>=nmax) return 0.0;
    return boost::math::laguerre(n, 1, 2.0*r/rscl) * exp(-r/rscl) / norm[n];
  } 

  //! Evaluate the the orthogonal Laguerre polynomial
  double eval(double& r, unsigned& n)
  {
    if (n>=nmax) return 0.0;
    return boost::math::laguerre(n, 1, 2.0*r/rscl) * exp(-r/rscl) / norm[n];
  } 

  //! Evaluate the the orthogonal Laguerre polynomial
  std::vector<double> eval(double& r)
  {
    std::vector<double> ret(nmax);

    // Initialization
    //
    double x = 2.0*r/rscl;
    ret[0] = boost::math::laguerre(0, 1, x);
    if (nmax>1) ret[1] = boost::math::laguerre(1, 1, x);

    // Recursion
    //
    for (int n=2; n<nmax; n++) 
      ret[n] = boost::math::laguerre_next(n-1, 1, x, ret[n-1], ret[n-2]);

    // Normalization
    //
    double w = exp(-r/rscl);
    for (int n=0; n<nmax; n++) ret[n] *= w / norm[n];

    return ret;
  } 

}; 


//! Coefficient file header
struct LaguCoefHeader
{
  const unsigned magic = 0x501acf;
  double time;
  double rscl;
  int nmax;
  int mnum;
};

class LaguCoefs
{
 private:

  double time, rscl, maccum;
  int nmax, mmax;

  using LaguPtr = std::shared_ptr<Laguerre>;

  std::map<int, LaguPtr> lagu;

public:

  //! Coefficient data
  std::map<int, std::array<std::vector<double>, 4>> cos_c, sin_c;

  //! Constructor
  LaguCoefs(double time, double rscl, int mmax, unsigned nmax) :
    time(time), rscl(rscl), mmax(mmax), nmax(nmax)
  {
    // Zero all accumulators
    //
    maccum = 0.0;
    for (int m=0; m<=mmax; m++) {
      lagu[m] = std::make_shared<Laguerre>(rscl, nmax);
      for (size_t k=0; k<4; k++) {
	cos_c[m][k].resize(nmax, 0);
	if (m) sin_c[m][k].resize(nmax, 0);
      }
    }
  }

  //! Add a particle contribution to coefficient
  void add(double mass, double R, double phi, double vr, double vt, double vz);

  //! MPI synchronize
  void synchronize()
  {
    MPI_Allreduce(MPI_IN_PLACE, &maccum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (int m=0; m<=mmax; m++) {
      for (size_t k=0; k<4; k++) {
	MPI_Allreduce(MPI_IN_PLACE, cos_c[m][k].data(), nmax, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	if (m)
	  MPI_Allreduce(MPI_IN_PLACE, sin_c[m][k].data(), nmax, MPI_DOUBLE,
			MPI_SUM, MPI_COMM_WORLD);
      }
      // END: k=0,...,4
    }
    // END: m loop
  }

  //! Write binary file
  void write(std::ostream& out);
};

typedef std::shared_ptr<LaguCoefs> LaguCoefPtr;

void LaguCoefs::write(std::ostream& out)
{
  LaguCoefHeader header;

  header.time   = time;
  header.rscl   = rscl;
  header.nmax   = nmax;
  header.mnum   = cos_c.size();
  
  out.write((const char *)&header, sizeof(LaguCoefHeader));

  for (auto d : cos_c) {
    out.write((const char *)&d.first, sizeof(int));
    for (int k=0; k<4; k++)
      out.write((const char *)d.second[k].data(), sizeof(double)*nmax);

    if (d.first) {
      for (int k=0; k<4; k++)
	out.write((const char *)sin_c[d.first][k].data(), sizeof(double)*nmax);
    }
  }
}

void
LaguCoefs::add(double mass, double R, double phi, double vr, double vt, double vz)
{
  // Add to grid
  maccum += mass;

  for (int m=0; m<=mmax; m++) {

    double cosm  = std::cos(phi*m), sinm = std::sin(phi*m);

    std::vector<double> val = lagu[m]->eval(R);
    
    for (unsigned int n=0; n<nmax; n++) {

      // Angular normalization and mass weighting
      //
      double fact  = mass * val[n] * 0.5*M_2_SQRTPI;
      if (m==0) fact *= M_SQRT1_2;
    
      cos_c[m][0][n] += fact*cosm;
      cos_c[m][1][n] += fact*vr*cosm;
      cos_c[m][2][n] += fact*vt*cosm;
      cos_c[m][3][n] += fact*vz*cosm;
      if (m) {
	sin_c[m][0][n] += fact*sinm;
	sin_c[m][1][n] += fact*vr*sinm;
	sin_c[m][2][n] += fact*vt*sinm;
	sin_c[m][3][n] += fact*vz*sinm;
      }
    }
  }
}

				// Globals for exputil library
				// Unused here
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

int
main(int ac, char **av)
{
  //===================
  // MPI preliminaries
  //===================

  MPI_Init(&ac, &av);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &proc_namelen);

  char *prog = av[0];
  bool verbose = false, finegrain= false;
  std::string cname, tname, new_dir, suffix, work_dir;
  int axis, nmax, comp, mmax, ibeg, iend;
  double rscl;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("finegrain",       "fine-grained progress report")
    ("append",          "append to existing output file")
    ("beg,i",	        po::value<int>(&ibeg)->default_value(0),
     "initial snapshot index")
    ("end,e",	        po::value<int>(&iend)->default_value(std::numeric_limits<int>::max()),
     "final snapshot index")
    ("mmax,M",	        po::value<int>(&mmax)->default_value(4),
     "maximum Fourier component in bin")
    ("rscale,a",        po::value<double>(&rscl)->default_value(0.01),
     "exponential disk scale")
    ("nmax,n",	        po::value<int>(&nmax)->default_value(8),
     "maximum Laguerre order")
    ("name,c",	        po::value<std::string>(&cname)->default_value("comp"),
     "component name")
    ("dir,d",           po::value<std::string>(&new_dir)->default_value("./"),
     "rewrite directory location for SPL files")
    ("work,w",          po::value<std::string>(&work_dir)->default_value("."),
     "working directory for output file")
    ("type,t",          po::value<std::string>(&tname)->default_value("OUT"),
     "PSP output type (OUT or SPL)")
    ("runtag,T",        po::value<std::string>(&runtag)->default_value("run0"),
     "Runtag id")
    ("suffix,s",        po::value<std::string>(&suffix)->default_value("ring_coefs"),
     "Output file suffix")
    ;

  po::variables_map vm;

  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0) {
      std::cout << desc << std::endl;
      std::cout << "Example: " << std::endl;
      std::cout << "\t" << av[0]
		<< " --runtag=run001" << std::endl;
    }
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  if (vm.count("finegrain")) {
    finegrain = true;
  }

  int n;
  for (n=ibeg; n<=iend; n++) {

    std::ostringstream fname;
    fname << tname << "." << runtag << "."
	  << std::setw(5) << std::setfill('0') << n;

    std::string file = fname.str();

    if (!FileExists(file)) {
      if (myid==0) {
	std::cerr << "Error opening file <" << file << "> for input"
		  << std::endl;
      }
      break;
    }
  }
  iend = n-1;
  if (myid==0) {
    std::cerr << "Assuming last file has index <" << iend << ">"
	      << std::endl;
  }

  auto iosflags = ios::out;
  if (vm.count("append")) iosflags |= ios::app;

  std::string outcoefs = work_dir + "/" + runtag + "." + suffix;
  std::ofstream out(outcoefs, iosflags);
  if (!out) {
    if (myid==0) {
      std::cerr << "Error opening file <" << outcoefs << "> for output"
		<< std::endl;
    }
    exit(-1);
  }

  std::shared_ptr<boost::progress_display> progress;
  if (myid==0 and not verbose and not finegrain) {
    progress = std::make_shared<boost::progress_display>(iend - ibeg + 1);
  }

  for (int n=ibeg; n<=iend; n++) {

    std::ostringstream fname;
    fname << tname << "." << runtag << "."
	  << std::setw(5) << std::setfill('0') << n;

    std::string file = fname.str();

    if (!FileExists(file)) {
      cerr << "Error opening file <" << file << "> for input\n";
      break;
    }

    if (myid==0) {
      if (verbose) cerr << "Using filename: " << file << endl;
    }

				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (file.find("SPL") != std::string::npos)
      psp = std::make_shared<PSPspl>(file, new_dir);
    else
      psp = std::make_shared<PSPout>(file);


				// Now write a summary
				// -------------------
    if (verbose and myid==0) {
      
      psp->PrintSummary(cerr);
    
      std::cerr << "\nPSP file <" << file << "> has time <" 
	   << psp->CurrentTime() << ">\n";
    }

				// Dump ascii for each component
				// -----------------------------
    std::vector<double> pos(3), vel(3);

				// Make the arrays
				// ---------------

    LaguCoefs lagu(psp->CurrentTime(), rscl, mmax, nmax);

    std::array<std::map<int, std::vector<float>>, 3> vel_c, vel_s;

    PSPstanza* stanza;
    SParticle* part;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      unsigned int icnt = 0;

      if (myid==0 and finegrain) {
	std::cout << "Using filename: " << file << std::endl;
	progress = std::make_shared<boost::progress_display>(stanza->comp.nbod/numprocs);
      }

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
	
	if (icnt++ % numprocs) continue;

	// Cylindrical radius
	//
	double R = 0.0;
	for (int k=0; k<2; k++) R += part->pos(k) * part->pos(k);
	R = sqrt(R);

	double mass = part->mass();
	
	// Make cylindrical velocity bins
	//
	double phi  = std::atan2(part->pos(1), part->pos(0));
	double cosp = std::cos(phi);
	double sinp = std::sin(phi);

	// uvec vr:  cos(phi), sin(phi)
	double vr = cosp*part->vel(0) + sinp*part->vel(1);

	// uvec vt: -sin(phi), cos(phi)
	double vt = -sinp*part->vel(0) + cosp*part->vel(1);
	
	// uvec vz
	double vz = part->vel(2);

	// Add to grid
	lagu.add(mass, R, phi, vr, vt, vz);

	if (myid==0 and finegrain) ++(*progress);
      }
      // END: particle loop
    }
    // END: stanza loop

    // Prepare for output
    //
    lagu.synchronize();
    if (myid==0) {
      lagu.write(out);
    }

    if (myid==0 and not verbose and not finegrain) {
      ++(*progress);
    }
  }
  // END: OUT file loop
  
  if (myid==0) std::cout << std::endl;

  MPI_Finalize();

  return 0;
}
