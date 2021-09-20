#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

// Boost stuff
//
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace po = boost::program_options;

#include <localmpi.H>
#include <SLGridMP2.H>
#include <gaussQ.H>

int main(int argc, char** argv)
{
  bool use_mpi, use_logr;
  double scale, rmin, rmax, rs, dfac;
  int numr, cmap, diverge, Lmax, nmax, knots;
  std::string filename, cachefile;

  //====================
  // Parse command line
  //====================

  po::options_description desc("Check the consistency a spherical SL basis\nAllowed options");
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("mpi",                 po::value<bool>(&use_mpi)->default_value(false),
     "using parallel computation")
    ("logr",                po::value<bool>(&use_logr)->default_value(false),
     "logarithmic spacing for orthogonality check")
    ("cmap",                po::value<int>(&cmap)->default_value(0),
     "coordinates in SphereSL: use mapped (1) or linear(0) coordinates")
    ("scale",               po::value<double>(&scale)->default_value(1.0),
     "scaling from real coordinates to table")
    ("Lmax",                po::value<int>(&Lmax)->default_value(2),
     "maximum number of angular harmonics in the expansion")
    ("nmax",                po::value<int>(&nmax)->default_value(10),
     "maximum number of radial harmonics in the expansion")
    ("numr",                po::value<int>(&numr)->default_value(1000),
     "radial knots for the SL grid")
    ("rmin",                po::value<double>(&rmin)->default_value(0.0001),
     "minimum radius for the SL grid")
    ("rmax",                po::value<double>(&rmax)->default_value(1.95),
     "maximum radius for the SL grid")
    ("rs",                  po::value<double>(&rs)->default_value(0.067),
     "cmap scale factor")
    ("diverge",             po::value<int>(&diverge)->default_value(0),
     "cusp divergence for spherical model")
    ("dfac",                po::value<double>(&dfac)->default_value(1.0),
     "cusp divergence exponent for spherical model")
    ("knots",               po::value<int>(&knots)->default_value(40),
     "Number of Legendre integration knots")
    ("filename",            po::value<string>(&filename)->default_value("SLGridSph.model"),
     "model file")
    ("cache",               po::value<string>(&cachefile)->default_value(".slgrid_sph_cache"),
     "cache file")
    ;

  //===================
  // Parse options
  //===================
  po::variables_map vm;
  
  // Parse command line for control and critical parameters
  //
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error on command line: "
			   << e.what() << std::endl;
    return -1;
  }
  
  //===================
  // MPI preliminaries 
  //===================
  if (use_mpi) {
    local_init_mpi(argc, argv);
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << desc << std::endl << std::endl;
    }
    if (use_mpi) MPI_Finalize();
    return 1;
  }

  if (use_mpi) {
    SLGridSph::mpi = 1;		// Turn on MPI
  } else {
    SLGridSph::mpi = 0;		// Turn off MPI
  }

				// Set model file
  SLGridSph::model_file_name = filename;
  SLGridSph::sph_cache_name = cachefile;

				// Generate Sturm-Liouville grid
  auto ortho = boost::make_shared<SLGridSph>(Lmax, nmax, numr, rmin, rmax, 
					     true, cmap, rs, 0, 1.0, true);
  //                                         ^               ^       ^
  //                                         |               |       |
  // Use cache file--------------------------+               |       |
  //                                                         |       |
  // Model extrapoltion--------------------------------------+       |
  //                                                                 |
  // Turn on diagnostic output in SL creation------------------------+

				// Slaves exit
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




