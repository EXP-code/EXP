#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <localmpi.H>
#include <SLGridMP2.H>
#include <Eigen/Eigen>
#include <gaussQ.H>
#include <cxxopts.H>

int main(int argc, char** argv)
{
  bool use_mpi, logr = false, ortho_check = false, dbg = false;
  double A, scale, rmin, rmax, L;
  int numr, cmap, diverge, mmax, kmax, nmax, knots, num;
  std::string filename, cachefile, model;

  // MPI preliminaries 
  //
  {
    // Use this hack to detect an OpenMPI environment (not portable,
    // but it's the best workaround I have).
    //
    std::string found(getenv("OMPI_COMM_WORLD_SIZE"));
    if (found.size()) {
      local_init_mpi(argc, argv);
      use_mpi = true;
    }
  }

  // Parse command line
  //
  cxxopts::Options options(argv[0], "Check the consistency a spherical SL basis");

  options.add_options()
    ("h,help", "Print this help message")
    ("ortho", "Compute orthogonality matrix")
    ("logr", "Plot output grid with logarithmic spacing")
    ("debug", "Print debugging output")
    ("cmap", "coordinates in SphereSL: use mapped (1) or linear(0) coordinates",
     cxxopts::value<int>(cmap)->default_value("0"))
    ("scale", "scaling from real coordinates to table",
     cxxopts::value<double>(scale)->default_value("1.0"))
    ("A,length", "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("L,thick", "pillbox size",
     cxxopts::value<double>(L)->default_value("1.0"))
    ("mmax", "maximum number of angular harmonics in the expansion",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("nmax", "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
    ("kmax", "maximum number of vertical wave numbers in the expansion",
     cxxopts::value<int>(kmax)->default_value("10"))
    ("numr", "radial knots for the SL grid",
     cxxopts::value<int>(numr)->default_value("1000"))
    ("r,rmin", "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
    ("R,rmax", "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("20.0"))
    ("num", "number of output grid points",
     cxxopts::value<int>(num)->default_value("1000"))
    ("knots", "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("40"))
    ("cache", "cache file",
     cxxopts::value<std::string>(cachefile)->default_value(".slgrid_sph_cache"))
    ("model", "SL model target type",
     cxxopts::value<std::string>(model)->default_value("expon"))
    ;

  
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

  // Orthogonality?
  //
  if (vm.count("ortho")) ortho_check = true;

  // Log spacing?
  //
  if (vm.count("logr")) logr = true;

  // Debugging output?
  //
  if (vm.count("debug")) dbg = true;

  if (use_mpi) {
    SLGridCyl::mpi = 1;		// Turn on MPI
  } else {
    SLGridCyl::mpi = 0;		// Turn off MPI
  }

  SLGridCyl::A = A;		// Set default scale length

  // Generate Sturm-Liouville grid
  //
  auto ortho = std::make_shared<SLGridCyl>(mmax, nmax, numr, kmax, rmin, rmax,
					   L, true, cmap, scale, model, dbg);
  //                                       ^  ^     ^     ^      ^      ^
  //                                       |  |     |     |      |      |
  // Slab height --------------------------+  |     |     |      |      |
  //                                          |     |     |      |      |
  // Use model cache -------------------------+     |     |      |      |
  //                                                |     |      |      |
  // Coordinate mapping type -----------------------+     |      |      |
  //                                                      |      |      |
  // Radial scale size -----------------------------------+      |      |
  //                                                             |      |
  // Target density-potential for SL creation--------------------+      |
  //                                                                    |
  // Turn on diagnostic output in SL creation---------------------------+

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

  cout << "M, N, K? ";
  int M, K, N;
  std::cin >> M;
  std::cin >> N;
  std::cin >> K;

  M = std::max<int>(M, 0);
  M = std::min<int>(M, mmax);

  N = std::max<int>(N, 0);
  N = std::min<int>(N, nmax-1);

  K = std::max<int>(K, 0);
  K = std::min<int>(K, kmax);

  double ximin = ortho->r_to_xi(rmin);
  double ximax = ortho->r_to_xi(rmax);

  // For the radial output grid
  //
  double x, r, lrmin, lrmax;
  
  // BEGIN: file header
  //
  out << "# M=" << M << " N=" << N << " K=" << K << std::endl;

  out << "# "
      << std::setw(13) << " x |"
      << std::setw(15) << " r |"
      << std::setw(15) << " P |"
      << std::setw(15) << " D |"
      << std::endl;

  out << "# "
      << std::setw(13) << " [1] |"
      << std::setw(15) << " [2] |"
      << std::setw(15) << " [3] |"
      << std::setw(15) << " [4] |"
      << std::endl;

  out << "# "
      << std::setw(13) << std::setfill('-') << "+"
      << std::setw(15) << std::setfill('-') << "+"
      << std::setw(15) << std::setfill('-') << "+"
      << std::setw(15) << std::setfill('-') << "+"
      << std::endl << std::setfill(' ');
  //
  // END: file header

  // Choosing the radial spacing
  //
  double delr;
  if (logr and rmin>0.0) {
    lrmin = log(rmin);
    lrmax = log(rmax);
    delr = (lrmax - lrmin)/(num - 1);
  } else {
    logr = false;
  }

  // Print potential density pairs 
  //
  for (int k=0; k<num; k++) {

    if (logr) {
      r = exp(lrmin + delr*k);
      x = ortho->r_to_xi(r);
    } else {
      x = ximin + (ximax - ximin)*k/(num-1);
      r = ortho->xi_to_r(x);
    }
	    
    out << std::setw(15) << x << std::setw(15) << r
	<< std::setw(15) << ortho->get_pot (r, M, N, K)
	<< std::setw(15) << ortho->get_dens(r, M, N, K)
	<< std::endl;
  }

  // Compute the inner product of the pairs
  //
  if (ortho_check) {

    std::ofstream out("cyltest.ortho");

    LegeQuad lw(knots);

    Eigen::MatrixXd orth(nmax, nmax);
    orth.setZero();

    for (int k=0; k<knots; k++) {
      double xx =  ximin + (ximax - ximin)*lw.knot(k);
      double rr = ortho->xi_to_r(xx);
      double fac = lw.weight(k) * rr / ortho->d_xi_to_r(xx) * (ximax - ximin);

      for (int j=0; j<nmax; j++) {
	for (int l=0; l<nmax; l++) {
	  orth(j, l) += fac *
	    ortho->get_pot (rr, M, j, K) *
	    ortho->get_dens(rr, M, l, K);
	  if (std::isnan(ortho->get_pot (rr, M, j, K))) {
	    std::cout << "pot R=" << rr << std::endl;
	  }
	  if (std::isnan(ortho->get_dens (rr, M, l, K))) {
	    std::cout << "dens R=" << rr << std::endl;
	  }
	}
      }
    }

    out << orth << std::endl;
  }

  if (use_mpi) MPI_Finalize();

  return 0;
}
