/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Cross validation analysis for SphereSL
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

				// Boost stuff

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include "Particle.h"
#include <PSP2.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>
#include <foarray.H>

#include <localmpi.H>

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

int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, RMAX, rscale, minSNR0;
  int NICE, LMAX, NMAX, NSNR, NPART;
  int beg, end, stride, init, knots, num;
  std::string modelf, dir("./"), cname, prefix;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Cross-validation analysis for spherical models" << std::endl
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
    ("minSNR",              po::value<double>(&minSNR0)->default_value(0.01),
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
  // Make SL expansion
  // ==================================================

  SphericalModelTable halo(modelf);
  SphereSL::mpi  = true;
  SphereSL::NUMR = 4000;
  SphereSL ortho(&halo, LMAX, NMAX, 1, rscale, true, NPART);

  auto sl = ortho.basis();

  // ==================================================
  // Compute and table Q values
  // ==================================================

  double ximin = sl->r_to_xi(RMIN);
  double ximax = sl->r_to_xi(RMAX);

  // =================
  // Begin radial loop
  // =================

  LegeQuad lw(knots);

  std::map< std::pair<int, int>, std::shared_ptr<std::vector<double>> > Q;

  for (int L=0; L<=LMAX; L++) {

    for (int n=1; n<=NMAX; n++) {

      auto qq = std::make_shared<std::vector<double>>(num);

      for (int i=0; i<num; i++) {
	double x = ximin + (ximax-ximin)*i/(num-1);
	double r = sl->xi_to_r(x);

	// Q1
	//
	double Q1 = 0.0;
	for (int k=0; k<knots; k++) {
	  double xx =  ximin + (x - ximin)*lw.knot(k);
	  double rr = sl->xi_to_r(xx);
	  Q1 += lw.weight(k) * sl->get_dens(xx, L, n, 0) * pow(rr/r, 1.0+L) * rr / sl->d_xi_to_r(xx);
	}
	Q1 *= (x - ximin)/(2.0*L+1.0);
    

	// Q2
	//
	double Q2 = 0.0;
	for (int k=0; k<knots; k++) {
	  double xx =  x + (ximax - x)*lw.knot(k);
	  double rr = sl->xi_to_r(xx);
	  Q2 += lw.weight(k) * sl->get_dens(xx, L, n, 0) * pow(r/rr, L) * rr / sl->d_xi_to_r(xx);
	}
	Q2 *= (ximax - x)/(2.0*L+1.0);
    
	(*qq)[i] = Q1 + Q2;
      }

      Q[std::pair<int, int>(L, n)] = qq;
    }
  }


  // ==================================================
  // Phase space output loop
  // ==================================================

  std::string file;

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
    ortho.reset_coefs();

    SParticle *p = psp->GetParticle();
    int icnt = 0;
    do {
      if (icnt++ % numprocs == myid)
	ortho.accumulate(p->pos(0), p->pos(1), p->pos(0), p->mass());
      p = psp->NextParticle();
    } while (p);
    
    if (myid==0) std::cout << "done" << endl;
    
    //------------------------------------------------------------ 
      
    if (myid==0) cout << "Making coefficients . . . " << flush;
    ortho.make_coefs();
    if (myid==0) std::cout << "done" << endl;
    if (myid==0) cout << "Making covariance . . . " << flush;
    ortho.make_covar(verbose);
    if (myid==0) std::cout << "done" << endl;

    //------------------------------------------------------------ 
    

    std::vector<double> term1(LMAX+1);
    std::vector<double> term2(LMAX+1), work2(LMAX+1);
    std::vector<double> term3(LMAX+1), work3(LMAX+1);
    
				// Default SNR limits
    double minSNR = minSNR0;
    double maxSNR = ortho.getMaxSNR();
				// Sanity check
    double dx = (ximax - ximin)/(num - 1);

    if (maxSNR < minSNR) minSNR = maxSNR / 100.0;
    
				// Space SNR logarithmically?
    if (LOG) {
      minSNR = log(minSNR);
      maxSNR = log(maxSNR);
    }

    double dSNR = (maxSNR - minSNR)/(NSNR - 1);

    if (myid==0) {
      std::cout << "minSNR=" << minSNR << " maxSNR=" << maxSNR << " dSNR=" << dSNR << std::endl;
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
      auto coefs = ortho.get_trimmed(snr, Hall);

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
	for (int L=0; L<=LMAX; L++) {
	  int lbeg = L*L;	// Offset into coefficient array
	  for (int M=0; M<=L; M++) {
	    for (int n=1; n<=NMAX; n++) {
	      if (M==0)
		term1[L] += coefs[lbeg][n]*coefs[lbeg][n];
	      else {
		int ll = lbeg + 2*(M-1) + 1;
		term1[L] +=
		  coefs[ll+0][n]*coefs[ll+0][n] +
		  coefs[ll+1][n]*coefs[ll+1][n];
	      }
	    }
	  }
	}
      }

      // Term 2 and Term 3
      //
      for (int L=0; L<=LMAX; L++) {

	int lbeg = L*L;	// Offset into coefficient array
	
	// Particle loop
	//
	p = psp->GetParticle();
	int icnt = 0;
	do {
	  if (icnt++ % numprocs == myid) {
	    double r = 0.0, costh = 0.0, phi = 0.0;
	    for (int k=0; k<3; k++) r += p->pos(k)*p->pos(k);
	    r = sqrt(r);
	    if (r>0.0) costh = p->pos(2)/r;
	    phi = atan2(p->pos(1), p->pos(0));
	    double mass = p->mass();
	    
	    double x = sl->r_to_xi(r);
	    x = std::max<double>(ximin, x);
	    x = std::min<double>(ximax, x);
	    
	    int indx = floor( (x - ximin)/dx );
	    if (indx<0) indx = 0;
	    if (indx>=num-1) indx = num-2;
	    
	    double A = (ximin + dx*(indx+1) - x)/dx;
	    double B = (x - ximin - dx*(indx+0))/dx;
	    
	    for (int M=0; M<=L; M++) {

	      double ylm = Ylm01(L, M) * plgndr(L, M, costh);

	      for (int n=1; n<=NMAX; n++) {

		std::pair<int, int> I(L, n);
		double Qval = A*(*Q[I])[indx] + B*(*Q[I])[indx+1];
		double potl = sl->get_pot(x, L, n, 0);
		  
		if (M==0) {
		  work2[L] += mass*coefs[lbeg+0][n]*ylm*potl;
		  work3[L] += mass*coefs[lbeg+0][n]*ylm*Qval;
		} else {
		  int ll = lbeg + 2*(M-1) + 1;
		  double fac = (coefs[ll][n]*cos(phi*M) + coefs[ll+1][n]*sin(phi*M))*ylm;
		  work2[L] += mass*potl * fac;
		  work3[L] += mass*Qval * fac;
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

