#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

#include <localmpi.H>
#include <SLGridMP2.H>
#include <gaussQ.H>
#include <cxxopts.H>

int main(int argc, char** argv)
{
  bool use_mpi, use_logr;
  double scale, rmin, rmax, rs, dfac;
  int numr, cmap, diverge, Lmax, nmax, knots;
  std::string filename, cachefile;

  //====================
  // Parse command line
  //====================

  cxxopts::Options options(argv[0], "Check the consistency a spherical SL basis");

  options.add_options()
   ("h,help", "Print this help message")
   ("mpi", "using parallel computation",
     cxxopts::value<bool>(use_mpi)->default_value("false"))
   ("logr", "logarithmic spacing for orthogonality check",
     cxxopts::value<bool>(use_logr)->default_value("false"))
   ("cmap", "coordinates in SphereSL: use mapped (1) or linear(0) coordinates",
     cxxopts::value<int>(cmap)->default_value("0"))
   ("scale", "scaling from real coordinates to table",
     cxxopts::value<double>(scale)->default_value("1.0"))
   ("Lmax", "maximum number of angular harmonics in the expansion",
     cxxopts::value<int>(Lmax)->default_value("2"))
   ("nmax", "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
   ("numr", "radial knots for the SL grid",
     cxxopts::value<int>(numr)->default_value("1000"))
   ("rmin", "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
   ("rmax", "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("1.95"))
   ("rs", "cmap scale factor",
     cxxopts::value<double>(rs)->default_value("0.067"))
   ("diverge", "cusp divergence for spherical model",
     cxxopts::value<int>(diverge)->default_value("0"))
   ("dfac", "cusp divergence exponent for spherical model",
     cxxopts::value<double>(dfac)->default_value("1.0"))
   ("knots", "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("40"))
   ("filename", "model file",
     cxxopts::value<string>(filename)->default_value("SLGridSph.model"))
   ("cache", "cache file",
     cxxopts::value<string>(cachefile)->default_value(".slgrid_sph_cache"))
    ;

  
  //===================
  // MPI preliminaries 
  //===================
  if (use_mpi) {
    local_init_mpi(argc, argv);
  }

  //===================
  // Parse options
  //===================

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    if (use_mpi) MPI_Finalize();
    return 2;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << options.help() << std::endl << std::endl;
    }
    if (use_mpi) MPI_Finalize();
    return 1;
  }

  if (use_mpi) {
    SLGridSph::mpi = 1;		// Turn on MPI
  } else {
    SLGridSph::mpi = 0;		// Turn off MPI
  }
				// Generate Sturm-Liouville grid
  auto ortho = std::make_shared<SLGridSph>(filename, Lmax, nmax, numr, rmin, rmax, 
					   true, cmap, rs, 0, 1.0, cachefile, true);
  //                                       ^               ^                  ^
  //                                       |               |                  |
  // Use cache file------------------------+               |                  |
  //                                                       |                  |
  // Model extrapoltion------------------------------------+                  |
  //                                                                          |
  // Turn on diagnostic output in SL creation---------------------------------+

				// Workers exit
  if (use_mpi && myid>0) {
    MPI_Finalize();
    exit(0);
  }

  std::cout << "Filename? ";
  std::cin  >> filename;
  std::ofstream out (filename.c_str());
  if (!out) {
    std::cout << "Can't open <" << filename << "> for output" << endl;
    exit(-1);
  }

  cout << "Number of points? ";
  int num;
  cin >> num;

  cout << "L, N? ";
  int L, N;
  cin >> L;
  cin >> N;

  N = std::max<int>(N, 0);
  N = std::min<int>(N, nmax);

  double ximin = ortho->r_to_xi(rmin);
  double ximax = ortho->r_to_xi(rmax);

  double x, r, lrmin, lrmax;
  
  if (use_logr) {
    if (rmin<1.0e-16) use_logr = false;
    else {
      lrmin = log(rmin);
      lrmax = log(rmax);
    }
  }

  // ==================
  // BEGIN: file header
  // ==================

  out << "# "
      << std::setw(13) << " x |"
      << std::setw(15) << " r |"
      << std::setw(15) << " Q1 |"
      << std::setw(15) << " Q2 |"
      << std::setw(15) << " Q |"
      << std::endl;

  out << "# "
      << std::setw(13) << " [1] |"
      << std::setw(15) << " [2] |"
      << std::setw(15) << " [3] |"
      << std::setw(15) << " [4] |"
      << std::setw(15) << " [5] |"
      << std::endl;

  out << "# "
      << std::setw(13) << std::setfill('-') << "+"
      << std::setw(15) << std::setfill('-') << "+"
      << std::setw(15) << std::setfill('-') << "+"
      << std::setw(15) << std::setfill('-') << "+"
      << std::setw(15) << std::setfill('-') << "+"
      << std::endl << std::setfill(' ');

  // ================
  // END: file header
  // ================

  // =================
  // Begin radial loop
  // =================

  LegeQuad lw(knots);
  double lr = 0.0;

  for (int i=0; i<num; i++) {

    if (use_logr) {
      lr = lrmin + (lrmax-lrmin)*i/(num-1);
      x = ortho->r_to_xi(exp(lr));
    } else {
      x = ximin + (ximax-ximin)*i/(num-1);
    }

    r = ortho->xi_to_r(x);
    out << setw(15) << x
	<< setw(15) << r;

    // Q1
    //
    double Q1 = 0.0;
    if (use_logr) {
      double x0 = ortho->r_to_xi(exp(lrmin));
      for (int k=0; k<knots; k++) {
	double rr = exp(lrmin + (lr-lrmin)*lw.knot(k));
	double xx = ortho->r_to_xi(rr);
	Q1 += lw.weight(k) * ortho->get_dens(xx, L, N, 0) * pow(rr/r, 1.0+L) * rr*rr;
      }
      Q1 *= (lr - lrmin)/(2.0*L+1.0);
    } else {
      for (int k=0; k<knots; k++) {
	double xx =  ximin + (x - ximin)*lw.knot(k);
	double rr = ortho->xi_to_r(xx);
	Q1 += lw.weight(k) * ortho->get_dens(xx, L, N, 0) * pow(rr/r, 1.0+L) * rr / ortho->d_xi_to_r(xx);
      }
      Q1 *= (x - ximin)/(2.0*L+1.0);
    }
    

    // Q2
    //
    double Q2 = 0.0;
    if (use_logr) {
      for (int k=0; k<knots; k++) {
	double rr = exp(lr + (lrmax-lr)*lw.knot(k));
	double xx = ortho->r_to_xi(rr);
	Q2 += lw.weight(k) * ortho->get_dens(xx, L, N, 0) * pow(r/rr, L)  * rr*rr;
      }
      Q2 *= (lrmax - lr)/(2.0*L+1.0);
    } else {
      for (int k=0; k<knots; k++) {
	double xx =  x + (ximax - x)*lw.knot(k);
	double rr = ortho->xi_to_r(xx);
	Q2 += lw.weight(k) * ortho->get_dens(xx, L, N, 0) * pow(r/rr, L) * rr / ortho->d_xi_to_r(xx);
      }
      Q2 *= (ximax - x)/(2.0*L+1.0);
    }
    
    out << std::setw(15) << Q1 << std::setw(15) << Q2 << std::setw(15) << Q1 + Q2 << endl;
  }

  if (use_mpi) MPI_Finalize();

  return 0;
}




