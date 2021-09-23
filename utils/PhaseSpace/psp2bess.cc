/*
  Separate a psp structure and make a kinematic Fourier coefficients
  series in Bessel functions

  MDWeinberg 05/20/20
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

#include <boost/program_options.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <Progress.H>

#include <mpi.h>

namespace po = boost::program_options;

//
// MPI variables
//
int numprocs, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];



//! Generate orthonormal Bessel functions of integral order
class Bess
{
private:

  unsigned int nroots;
  double order;
  std::vector<double> roots, norm;

public:

  //! Constructor: set the order and number of radial functions
  Bess(double order, unsigned int nroots) : order(order), nroots(nroots)
  {
    boost::math::cyl_bessel_j_zero(order, 1, nroots, std::back_inserter(roots));
    norm.resize(nroots);
    for (unsigned int m=0; m<nroots; m++) {
      double val = boost::math::cyl_bessel_j(order+1.0, roots[m]);
      norm[m] = sqrt(0.5*val*val);
    }
  }
  
  //! Get the norm for radial order m
  double getNorm(int m)
  {
    if (m>=nroots) return 0.0;
    else           return norm[m];
  }
  
  //! Evaluate the Bessel for x in [0, 1] for radial order m
  double operator()(double& x, const unsigned& m)
  {
    if (m>=nroots) return 0.0;
    return boost::math::cyl_bessel_j<double, double>(order, x*roots[m]) / norm[m];
  } 

  //! Evaluate the Bessel for x in [0, 1] for radial order m
  double eval(double& x, unsigned& m)
  {
    if (m>=nroots) return 0.0;
    return boost::math::cyl_bessel_j<double, double>(order, x*roots[m]) / norm[m];
  } 

}; 


//! Coefficient file header
struct BessCoefHeader
{
  const unsigned magic = 0x501ace;
  double time;
  double rmax;
  int nmax;
  int mnum;
};

class BessCoefs
{
 private:

  double time, rmax, maccum;
  int nmax, mmax;

  using BessPtr = std::shared_ptr<Bess>;

  std::map<int, BessPtr> bess;

public:

  //! Coefficient data
  std::map<int, std::array<std::vector<double>, 4>> cos_c, sin_c;

  //! Constructor
  BessCoefs(double time, double rmax, int mmax, unsigned nmax) :
    time(time), rmax(rmax), mmax(mmax), nmax(nmax)
  {
    // Zero all accumulators
    //
    maccum = 0.0;
    for (int m=0; m<=mmax; m++) {
      bess[m] = std::make_shared<Bess>(static_cast<double>(m), nmax);
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

typedef std::shared_ptr<BessCoefs> BessCoefPtr;

void BessCoefs::write(std::ostream& out)
{
  BessCoefHeader header;

  header.time   = time;
  header.rmax   = rmax;
  header.nmax   = nmax;
  header.mnum   = cos_c.size();
  
  out.write((const char *)&header, sizeof(BessCoefHeader));

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
BessCoefs::add(double mass, double R, double phi, double vr, double vt, double vz)
{
  if (R>rmax) return;

  // Add to grid
  maccum += mass;

  for (int m=0; m<=mmax; m++) {
    double cosm  = std::cos(phi*m), sinm = std::sin(phi*m);

    for (unsigned int n=0; n<nmax; n++) {

      double x     = R/rmax;
      double value = bess[m]->eval(x, n)/rmax;
      //                                 ^
      //                                 |
      // Scale factor in normalization---+

      // Angular normalization and mass weighting
      //
      double fact  = mass * value * 0.5*M_2_SQRTPI;
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
  std::string cname, tname, new_dir, suffix, work_dir, runtag;
  int axis, nmax, comp, mmax, ibeg, iend;
  double rmax;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("finegrain",       "fine-grained progress report")
    ("beg,i",	        po::value<int>(&ibeg)->default_value(0),
     "initial snapshot index")
    ("end,e",	        po::value<int>(&iend)->default_value(std::numeric_limits<int>::max()),
     "final snapshot index")
    ("mmax,M",	        po::value<int>(&mmax)->default_value(4),
     "maximum Fourier component in bin")
    ("rmax,R",	        po::value<double>(&rmax)->default_value(0.04),
     "maximum radius")
    ("nmax,n",	        po::value<int>(&nmax)->default_value(8),
     "maximum Bessel order")
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

  std::string outcoefs = work_dir + "/" + runtag + "." + suffix;
  std::ofstream out(outcoefs, ios::app | ios::out);
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

    BessCoefs bess(psp->CurrentTime(), rmax, mmax, nmax);

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
	bess.add(mass, R, phi, vr, vt, vz);

	if (myid==0 and finegrain) ++(*progress);
      }
      // END: particle loop
    }
    // END: stanza loop

    // Prepare for output
    //
    bess.synchronize();
    if (myid==0) {
      bess.write(out);
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
