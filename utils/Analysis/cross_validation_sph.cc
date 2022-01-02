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

#include <global.H>
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
  std::string modelf, dir("./"), cname, prefix, fileType, filePrefix;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  cxxopts::Options options("pspxvalH",  "Cross-validation analysis for spherical models");

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
     cxxopts::value<string>(prefix)->default_value("crossval"))
    ("runtag", "Phase space file",
     cxxopts::value<string>(runtag)->default_value("run1"))
    ("outdir", "Output directory path",
     cxxopts::value<string>(outdir)->default_value("."))
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
	    sl->get_dens(xx, L, n1, 0) *sl->get_pot(xx, L, n2, 0) * rr * rr;
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

    auto file1 = PR::ParticleReader::fileNameCreator
      (fileType, ipsp, myid, dir, runtag);
    
    PR::PRptr reader = PR::ParticleReader::createReader
      (fileType, file1, true);

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
	ortho.accumulate(p->pos[0], p->pos[1], p->pos[0], p->mass);
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
    std::vector<double> term2(LMAX+1), work2(LMAX+1);
    std::vector<double> term3(LMAX+1), work3(LMAX+1);
    
				// Sanity check
    double dx = (ximax - ximin)/(num - 1);

    if (NCUT) {

      double term4tot = 0.0;

      auto coefs = ortho.retrieve_coefs();

      for (int ncut=0; ncut<NMAX; ncut++) {

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
	      for (int n1=0; n1<ncut; n1++) {
		for (int n2=0; n2<ncut; n2++) {
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
	for (int L=0; L<=LMAX; L++) {
	  
	  int lbeg = L*L;	// Offset into coefficient array
	  
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
	      
	      double x = sl->r_to_xi(r);
	      x = std::max<double>(ximin, x);
	      x = std::min<double>(ximax, x);
	      
	      int indx = floor( (x - ximin)/dx );
	      if (indx<0) indx = 0;
	      if (indx>=num-1) indx = num-2;
	      
	      double A = (ximin + dx*(indx+1) - x)/dx;
	      double B = (x - ximin - dx*(indx+0))/dx;
	      
	      for (int M=0; M<=L; M++) {

		double ylm = Ylm_fac(L, M) * plgndr(L, M, costh);

		for (int n=0; n<ncut; n++) {

		  std::pair<int, int> I(L, n);
		  double Qval = A*(*Q[I])[indx] + B*(*Q[I])[indx+1];
		  double potl = sl->get_pot(x, L, n, 0);
		
		  if (M==0) {
		    if (r<RMAX)
		      work2[L] += mass*coefs(lbeg+0, n)*ylm*potl;
		    work3[L] += -mass*coefs(lbeg+0, n)*ylm*Qval;
		  } else {
		    int ll = lbeg + 2*(M-1) + 1;
		    double fac = (coefs(ll, n)*cos(phi*M) + coefs(ll+1, n)*sin(phi*M))*ylm;
		    if (r<RMAX) work2[L] += mass*potl * fac;
		    work3[L] += -mass*Qval * fac;
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
	  
	  out << std::setw( 5) << ipsp
	      << std::setw( 5) << ncut;
	
	  double term1tot = std::accumulate(term1.begin(), term1.end(), 0.0) / (4.0*M_PI);
	  double term2tot = std::accumulate(term2.begin(), term2.end(), 0.0);
	  double term3tot = std::accumulate(term3.begin(), term3.end(), 0.0);
	  
	  if (ncut==0) term4tot = term1tot;
	  
	  out << std::setw(18) << term1tot
	      << std::setw(18) << term2tot
	      << std::setw(18) << term3tot
	      << std::setw(18) << term1tot - term2tot - term3tot + term4tot
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
	for (int L=0; L<=LMAX; L++) {
	  
	  int lbeg = L*L;	// Offset into coefficient array
	
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
	    
	      double x = sl->r_to_xi(r);
	      x = std::max<double>(ximin, x);
	      x = std::min<double>(ximax, x);
	    
	      int indx = floor( (x - ximin)/dx );
	      if (indx<0) indx = 0;
	      if (indx>=num-1) indx = num-2;
	    
	      double A = (ximin + dx*(indx+1) - x)/dx;
	      double B = (x - ximin - dx*(indx+0))/dx;
	    
	      for (int M=0; M<=L; M++) {

		double ylm = Ylm_fac(L, M) * plgndr(L, M, costh);

		for (int n=1; n<=NMAX; n++) {

		  std::pair<int, int> I(L, n);
		  double Qval = A*(*Q[I])[indx] + B*(*Q[I])[indx+1];
		  double potl = sl->get_pot(x, L, n, 0);
		  
		  if (M==0) {
		    if (r<RMAX)
		      work2[L] += mass*coefs(lbeg+0, n)*ylm*potl;
		    work3[L] += mass*coefs(lbeg+0, n)*ylm*Qval;
		  } else {
		    int ll = lbeg + 2*(M-1) + 1;
		    double fac = (coefs(ll, n)*cos(phi*M) + coefs(ll+1, n)*sin(phi*M))*ylm;
		    if (r<RMAX)
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

