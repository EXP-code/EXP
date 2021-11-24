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


                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>
#include <foarray.H>

#include <localmpi.H>
#include <cxxopts.H>

// Globals
//

extern double plgndr(int, int, double);

double Ylm_fac(int ll, int mm)
{
  mm = abs(mm);
  return sqrt( (2.0*ll+1)/(4.0*M_PI*M_PI) ) *
    exp(0.5*(lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm)));
}

int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, RMAX, rscale, minSNR0, Hexp;
  int NICE, LMAX, NMAX, NSNR, NPART;
  int beg, end, stride, init, knots, num;
  std::string modelf, dir("./"), cname, prefix, fileType, filePrefix, runtag;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Cross-validation analysis for spherical models" << std::endl
       << std::string(60, '-') << std::endl << std::endl
       << "Allowed options";
  
  cxxopts::Options options("pspxval2", sout.str());

  options.add_options()
    ("h,help", "Print this help message")
    ("v,verbose", "Verbose and diagnostic output for covariance computation")
    ("NCUT", "trim coefficient by order rather than SNR")
    ("LOG", "log scaling for SNR")
    ("Hall", "use Hall smoothing for SNR trim")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("P,prefix", "prefix for phase-space files",
     cxxopts::value<std::string>(filePrefix)->default_value("OUT"))
    ("NICE", "system priority",
     cxxopts::value<int>(NICE)->default_value("0"))
    ("RMIN", "minimum radius for Q table",
     cxxopts::value<double>(RMIN)->default_value("0.0"))
    ("RSCALE", "coordinate mapping scale factor",
     cxxopts::value<double>(rscale)->default_value("0.067"))
    ("RMAX", "maximum radius for Q table",
     cxxopts::value<double>(RMAX)->default_value("2.0"))
    ("LMAX", "Maximum harmonic order for spherical expansion",
     cxxopts::value<int>(LMAX)->default_value("4"))
    ("NMAX", "Maximum radial order for spherical expansion",
     cxxopts::value<int>(NMAX)->default_value("12"))
    ("NPART", "Jackknife partition number for testing (0 means off, use standard eval)",
     cxxopts::value<int>(NPART)->default_value("0"))
    ("N,NSNR", "Number of SNR evaluations",
     cxxopts::value<int>(NSNR)->default_value("20"))
    ("minSNR", "minimum SNR value for loop output",
     cxxopts::value<double>(minSNR0)->default_value("0.01"))
    ("Hexp", "default Hall smoothing exponent",
     cxxopts::value<double>(Hexp)->default_value("1.0"))
    ("prefix", "Filename prefix",
     cxxopts::value<std::string>(prefix)->default_value("crossval"))
    ("runtag", "Phase space file",
     cxxopts::value<string>(runtag)->default_value("run1"))
    ("modelfile", "Halo model file",
     cxxopts::value<string>(modelf)->default_value("SLGridSph.model"))
    ("init", "fiducial PSP index",
     cxxopts::value<int>(init)->default_value("0"))
    ("beg", "initial PSP index",
     cxxopts::value<int>(beg)->default_value("0"))
    ("end", "final PSP index",
     cxxopts::value<int>(end)->default_value("99999"))
    ("stride", "PSP index stride",
     cxxopts::value<int>(stride)->default_value("1"))
    ("num", "Number of entries in Q table",
     cxxopts::value<int>(num)->default_value("10000"))
    ("knots", "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("40"))
    ("compname", "train on Component (default=stars)",
     cxxopts::value<std::string>(cname)->default_value("stars"))
    ("d,dir", "directory for SPL files",
     cxxopts::value<std::string>(dir))
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
  if (vm.count("LOG")) LOG = true;

  bool Hall = false;
  if (vm.count("Hall")) Hall = true;

  bool verbose = false;
  if (vm.count("verbose")) verbose = true;

  bool NCUT = false;
  if (vm.count("NCUT")) NCUT = true;

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
  // Make SL expansion
  // ==================================================

  auto halo = std::make_shared<SphericalModelTable>(modelf);

  SphereSL::mpi  = true;
  SphereSL::NUMR = 4000;
  SphereSL::HEXP = Hexp;

  SphereSL ortho(halo, LMAX, NMAX, 1, rscale, true, NPART);

  auto sl = ortho.basis();

  // ==================================================
  // Compute and table M and Q values
  // ==================================================

  double ximin = sl->r_to_xi(RMIN);
  double ximax = sl->r_to_xi(RMAX);

  // =================
  // Begin radial loop
  // =================

  LegeQuad lw(knots);

  std::map<int, Eigen::MatrixXd> O;
  for (int L=0; L<=LMAX; L++) {
    O[L] = Eigen::MatrixXd::Zero(NMAX, NMAX);
    
    for (int n1=1; n1<=NMAX; n1++) {
      for (int n2=1; n2<=NMAX; n2++) {
	for (int k=0; k<knots; k++) {
	  double xx =  ximin + (ximax - ximin)*lw.knot(k);
	  double rr = sl->xi_to_r(xx);
	  O[L](n1-1, n2-1) +=
	    lw.weight(k) * (ximax - ximin) / sl->d_xi_to_r(xx) *
	    sl->get_dens(xx, L, n1, 0) *sl->get_dens(xx, L, n2, 0) * rr * rr;
	}
      }
    }
  }

  // Print orthogonality matrix
  //
  if (false and myid==0) {
    std::cout << std::string(60, '-') << std::endl;
    for (auto v : O) {
      std::cout << std::setw(4) << v.first << std::endl
		<< v.second << std::endl;
      std::cout << std::string(60, '-') << std::endl;
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
    auto file1 = ParticleReader::fileNameCreator
      (fileType, ipsp, myid, dir, runtag);

    std::ifstream in(file1);
    if (!in) {
      std::cerr << "Error opening <" << file << ">" << endl;
      iok = 0;
    }
    
    if (iok==0) break;

    // ==================================================
    // Open PSP file
    // ==================================================

    PRptr reader = ParticleReader::createReader(fileType, file1, true);

    double tnow = reader->CurrentTime();
    if (myid==0) std::cout << "Beginning partition [time=" << tnow
			   << ", index=" << ipsp << "] . . . "  << flush;
    
    reader->SelectType(cname);

    //------------------------------------------------------------ 

    if (myid==0) std::cout << std::endl
			   << "Accumulating particle positions . . . "
			   << std::flush;
    ortho.reset_coefs();

    auto p = reader->firstParticle();
    int icnt = 0;
    do {
      if (icnt++ % numprocs == myid)
	ortho.accumulate(p->pos[0], p->pos[1], p->pos[2], p->mass);
      p = reader->nextParticle();
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
    double term2, work2, pi42 = 0.0625/M_PI/M_PI;
    
				// Sanity check
    double dx = (ximax - ximin)/(num - 1);

    if (NCUT) {

      double term4tot = 0.0;

      for (int ncut=1; ncut<=NMAX; ncut++) {

	// Zero out the accumulators
	//
	std::fill(term1.begin(), term1.end(), 0.0);
	term2 = work2 = 0.0;
	
	auto coefs = ortho.retrieve_coefs();
	for (int l=0; l<=coefs.rows(); l++) {
	  for (int n=ncut+1; n<NMAX; n++) coefs(l, n) = 0.0;
	}

	if (myid==0) {		// Only root process needs this one
	  
	  // Term 1
	  //
	  double t1 = 0.0;
	  for (int L=0; L<=LMAX; L++) {
	    int lbeg = L*L;	// Offset into coefficient array
	    for (int M=0; M<=L; M++) {
	      for (int n1=1; n1<=ncut; n1++) {
		for (int n2=1; n2<=ncut; n2++) {
		  if (M==0)
		    term1[L] +=
		      coefs(lbeg, n1)*O[L](n1-1, n2-1)*coefs(lbeg, n2);
		  else {
		    int ll = lbeg + 2*(M-1) + 1;
		    term1[L] +=
		      coefs(ll+0, n1)*O[L](n1-1, n2-1)*coefs(ll+0, n2) +
		      coefs(ll+1, n1)*O[L](n1-1, n2-1)*coefs(ll+1, n2);
		  }
		}
	      }
	    }
	  }
	}
	
	// Term 2 and Term 3
	//
	  
	// Particle loop
	//
	p = reader->firstParticle();
	int icnt = 0;
	do {
	  if (icnt++ % numprocs == myid) {
	    double r = 0.0, costh = 0.0, phi = 0.0;
	    for (int k=0; k<3; k++) r += p->pos[k]*p->pos[k];
	    r = sqrt(r);
	    if (r>0.0) costh = p->pos[2]/r;
	    phi = atan2(p->pos[1], p->pos[0]);
	    double mass = p->mass;
	    
	    double d0, d, p0, p;
	    ortho.dens_pot_eval(r, costh, phi, d0, d, p0, p);
	    
	    work2 += mass*(d0 + d);
	  }
	  // END: add particle data
	    
	  // Queue up next particle
	  //
	  p = reader->nextParticle();
	} while (p);
	//
	// END: particle loop

	
	MPI_Reduce(&work2, &term2, 1, MPI_DOUBLE,
		   MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (myid==0) {
	  
	  out << std::setw( 5) << ipsp
	      << std::setw( 5) << ncut;
	  
	  double term1tot = std::accumulate(term1.begin(), term1.end(), 0.0) * pi42;
	  
	  if (ncut==1) term4tot = term1tot;
	  
	  out << std::setw(18) << term1tot
	      << std::setw(18) << term2
	      << std::setw(18) << term1tot - 2.0*term2 + term4tot
	      << std::endl;
	}
	// Root process
	
	if (myid==0) std::cout << "done" << endl;
	
      }
      // NCUT loop

    } else {

      double minSNR = minSNR0;
      double maxSNR = ortho.getMaxSNR();

      if (maxSNR < minSNR) minSNR = maxSNR / 100.0;
      
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
	auto coefs = ortho.get_trimmed(snr, ortho.getMass(), Hall);

	// 

	ortho.install_coefs(coefs);
	
	// Zero out the accumulators
	//
	std::fill(term1.begin(), term1.end(), 0.0);
	term2 = work2 = 0.0;
	
	if (myid==0) {		// Only root process needs this one
	  
	  // Term 1
	  //
	  double t1 = 0.0;
	  for (int L=0; L<=LMAX; L++) {
	    int lbeg = L*L;	// Offset into coefficient array
	    for (int M=0; M<=L; M++) {
	      for (int n1=1; n1<=NMAX; n1++) {
		for (int n2=1; n2<=NMAX; n2++) {
		  if (M==0)
		    term1[L] +=
		      coefs(lbeg, n1)*O[L](n1-1, n2-1)*coefs(lbeg, n2);
		  else {
		    int ll = lbeg + 2*(M-1) + 1;
		    term1[L] +=
		      coefs(ll+0, n1)*O[L](n1-1, n2-1)*coefs(ll+0, n2) +
		      coefs(ll+1, n1)*O[L](n1-1, n2-1)*coefs(ll+1, n2);
		  }
		}
	      }
	    }
	  }
	}
	
	// Term 2 and Term 3
	//
	
	// Particle loop
	//
	p = reader->firstParticle();
	int icnt = 0;
	do {
	  if (icnt++ % numprocs == myid) {
	    double r = 0.0, costh = 0.0, phi = 0.0;
	    for (int k=0; k<3; k++) r += p->pos[k]*p->pos[k];
	    r = sqrt(r);
	    if (r>0.0) costh = p->pos[2]/r;
	    phi = atan2(p->pos[1], p->pos[0]);
	    double mass = p->mass;
	    
	    double d0, d, p0, p;
	    ortho.dens_pot_eval(r, costh, phi, d0, d, p0, p);
	    
	    work2 += mass*(d0 + d);
	  }
	  // END: add particle data
	    
	  // Queue up next particle
	  //
	  p = reader->nextParticle();
	} while (p);
	//
	// END: particle loop
	
	MPI_Reduce(&work2, &term2, 1, MPI_DOUBLE,
		   MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (myid==0) {
	  
	  out << std::setw( 5) << ipsp
	      << std::setw(18) << snr;
	  
	  double term1tot = std::accumulate(term1.begin(), term1.end(), 0.0) * pi42;
	  
	  if (nsnr==0) term4tot = term1tot;
	  
	  out << std::setw(18) << term1tot
	      << std::setw(18) << term2
	      << std::setw(18) << term1tot - 2.0*term2 + term4tot
	      << std::endl;
	}
	// Root process

	if (myid==0) std::cout << "done" << endl;
      }
      // END: snr loop
    }
    // Method switch
      
    if (myid==0) std::cout << "done" << endl;

    // Blank line between stanzas
    //
    if (myid==0) out << std::endl;
  }
  // END: dump loop


  MPI_Finalize();

  return 0;
}

