/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Cross validation analysis for sphere
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
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

using namespace std;

				// Boost stuff

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
#include <SphereSL.H>
#include <localmpi.h>
#include <foarray.H>

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


int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, RMAX;
  int NICE, LMAX, NMAX, L1, L2;
  int beg, end, stride, init, knots, num;
  std::string MODFILE, INDEX, dir("./"), cname, outfile;

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
    ("OUT",
     "assume original, single binary PSP files as input")
    ("SPL",
     "assume new split binary PSP files as input")
    ("NICE",                po::value<int>(&NICE)->default_value(0),
     "system priority")
    ("RMIN",                po::value<double>(&RMIN)->default_value(0.0),
     "minimum radius for output")
    ("RMAX",                po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("LMAX",                po::value<int>(&LMAX)->default_value(4),
     "Maximum harmonic order for spherical expansion")
    ("NMAX",                po::value<int>(&NMAX)->default_value(12),
     "Maximum radial order for spherical expansion")
    ("L1",                  po::value<int>(&L1)->default_value(0),
     "minimum l harmonic")
    ("L2",                  po::value<int>(&L2)->default_value(100),
     "maximum l harmonic")
    ("OUTFILE",             po::value<string>(&outfile)->default_value("crossval.dat"),
     "Filename prefix")
    ("runtag",              po::value<string>(&runtag)->default_value("run1"),
     "Phase space file")
    ("outdir",              po::value<string>(&outdir)->default_value("."),
     "Output directory path")
    ("MODFILE",             po::value<string>(&MODFILE)->default_value("SLGridSph.model"),
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
    std::cout << std::endl << desc << std::endl;
    return 0;
  }

  bool SPL = false;
  if (vm.count("SPL")) SPL = true;
  if (vm.count("OUT")) SPL = false;

  // ==================================================
  // Nice process
  // ==================================================

  if (NICE>0)
    setpriority(PRIO_PROCESS, 0, NICE);

  std::ofstream out(outfile);
  if (not out) {
    std::cerr << "Error opening output file <" << outfile << ">" << std::endl;
    exit(-2);
  }

  // ==================================================
  // Make SL expansion
  // ==================================================

  SphericalModelTable halo(MODFILE);
  SphereSL::mpi = true;
  SphereSL::NUMR = 4000;
  SphereSL ortho(&halo, LMAX, NMAX, 1);

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

  std::map< std::pair<int, int>, std::vector<double> > Q;

  for (int L=L1; L<=L2; L++) {
    for (int n=1; n<=NMAX; n++) {

      std::vector<double> qq(NMAX);

      for (int i=0; i<num; i++) {
	double x = ximin + (ximax-ximin)*i/(num-1);
	double r = sl->xi_to_r(x);

	// Q1
	//
	double Q1 = 0.0;
	for (int k=1; k<=knots; k++) {
	  double xx =  ximin + (x - ximin)*lw.knot(k);
	  double rr = sl->xi_to_r(xx);
	  Q1 += lw.weight(k) * sl->get_dens(xx, L, n, 0) * pow(rr/r, 1.0+L) * rr / sl->d_xi_to_r(xx);
	}
	Q1 *= (x - ximin)/(2.0*L+1.0);
    

	// Q2
	//
	double Q2 = 0.0;
	for (int k=1; k<=knots; k++) {
	  double xx =  x + (ximax - x)*lw.knot(k);
	  double rr = sl->xi_to_r(xx);
	  Q2 += lw.weight(k) * sl->get_dens(xx, L, n, 0) * pow(r/rr, L) * rr / sl->d_xi_to_r(xx);
	}
	Q2 *= (ximax - x)/(2.0*L+1.0);
    
	qq[i] = Q1 + Q2;
      }

      Q[std::pair<int, int>(L, n)] = qq;
    }
  }


  // ==================================================
  // Phase space output loop
  // ==================================================

  std::string file;

  for (int indx=beg; indx<=end; indx+=stride) {

    // ==================================================
    // PSP input stream
    // ==================================================

    int iok = 1;
    std::ostringstream s0, s1;

    if (myid==0) {
      s1.str("");		// Clear stringstream
      if (SPL) s1 << "SPL.";
      else     s1 << "OUT.";
      s1 << runtag << "."<< std::setw(5) << std::setfill('0') << indx;
      
				// Check for existence of next file
      file = dir + s1.str();
      std::ifstream in(file);
      if (!in) {
	std::cerr << "Error opening <" << file << ">" << endl;
	iok = 0;
      }
    }
    
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (iok==0) break;

    // ==================================================
    // Open PSP file
    // ==================================================
    PSPptr psp;

    if (myid==0) {

      if (SPL) psp = std::make_shared<PSPspl>(s1.str(), dir, true);
      else     psp = std::make_shared<PSPout>(file, true);

      tnow = psp->CurrentTime();
      cout << "Beginning partition [time=" << tnow
	   << ", index=" << indx << "] . . . "  << flush;
    }
    
    if (not psp->GetNamed(cname)) {
      std::cout << "Error finding component named <" << cname << ">" << std::endl;
      psp->PrintSummary(std::cout);
      exit(-1);
    }
    
    //------------------------------------------------------------ 

    if (myid==0) cout << "Accumulating particle positions . . . " << flush;

    ortho.reset_coefs();

    SParticle *p = psp->GetParticle();
    do {
      ortho.accumulate(p->pos(0), p->pos(1), p->pos(0), p->mass());
      p = psp->NextParticle();
    } while (p);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
  
    //------------------------------------------------------------ 

    if (myid==0) cout << "Making coefficients . . . " << flush;
    ortho.make_coefs();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    double time = 0.0;
    if (myid==0) {
      double term1 = 0.0;
      double term2 = 0.0;
      double term3 = 0.0;

      auto coefs = ortho.retrieve_coefs();

      for (int n=1; n<=NMAX; n++) {

	// Term 1
	//
	double t1 = 0.0;
	for (int L=L1; L<=L2; L++) {
	  int lbeg = L1*L1;
	  for (int M=0; M<=L; M++) {
	    if (M==0)
	      term1 += coefs[lbeg][n]*coefs[L][n];
	    else {
	      int ll = lbeg + 2*M + 1;
	      term1 +=
		coefs[lbeg+0][n]*coefs[lbeg+0][n] +
		coefs[lbeg+1][n]*coefs[lbeg+1][n];
	    }
	  }
	}
	  

	// Term 2
	//
	for (int L=L1; L<=L2; L++) {
	  int lbeg = L1*L1;

	  p = psp->GetParticle();
	  do {
	    double r = 0.0, costh = 0.0, phi = 0.0;
	    for (int k=0; k<3; k++) r += p->pos(k)*p->pos(k);
	    r = sqrt(r);
	    if (r>0.0) costh = p->pos(2)/r;
	    phi = atan2(p->pos(1), p->pos(0));
	    double mass = p->mass();

	    double x = sl->r_to_xi(r);
	    x = std::max<double>(ximin, x);
	    x = std::min<double>(ximax, x);
	    
	    double dens = sl->get_dens(x, L, n, 0);

	    int indx = (x - ximin)/(ximax - ximin)*num;
	    if (indx<0) indx = 1;
	    if (indx>=NMAX-1) indx = NMAX-2;
	    
	    double A = (ximin + (ximax - ximin)*(indx+1)/num - x)/(ximax - ximin);
	    double B = (x - ximin - (ximax - ximin)*(indx+0)/num)/(ximax - ximin);
	    std::pair<int, int> I(L, n);
	    double Qval = A*Q[I][indx] + B*Q[I][indx+1];


	    for (int M=0; M<=L; M++) {
	      double ylm = Ylm01(L, M);
	      if (M==0) {
		term2 += mass*coefs[lbeg+0][n]*ylm*dens;
		term3 += mass*coefs[lbeg+0][n]*ylm*Qval;
	      } else {
		int ll = lbeg + 2*M + 1;
		double fac = (coefs[ll][n]*cos(phi*M) + coefs[ll+1][n]*sin(phi*M))*ylm;
		
		term2 += mass*dens * fac;
		term3 += mass*Qval * fac;
	      }
	    }
	    p = psp->NextParticle();
	  } while (p);
	}

	out << std::setw( 5) << beg
	    << std::setw( 5) << n
	    << std::setw(18) << term1
	    << std::setw(18) << term2
	    << std::setw(18) << term3
	    << std::endl;
      }

    }
    MPI_Barrier(MPI_COMM_WORLD);

    //------------------------------------------------------------ 

  } // Dump loop

  MPI_Finalize();

  return 0;
}

