#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

#include <Eigen/Eigen>

#include <libvars.H>
#include <localmpi.H>
#include <SLGridMP2.H>
#include <gaussQ.H>
#include <cxxopts.H>

int main(int argc, char** argv)
{
  bool use_mpi, use_logr;
  double rmin, rmax, rs, dfac;
  int numr, cmap, diverge, Lmax, nmax;
  string filename, cachefile;

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
   ("Lmax", "maximum number of angular harmonics in the expansion",
     cxxopts::value<int>(Lmax)->default_value("2"))
   ("nmax", "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
   ("numr", "radial knots for the SL grid",
     cxxopts::value<int>(numr)->default_value("1000"))
   ("rmin", "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("-1.0"))
   ("rmax", "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("-1.0"))
   ("rs", "cmap scale factor",
     cxxopts::value<double>(rs)->default_value("0.067"))
   ("diverge", "cusp divergence for spherical model",
     cxxopts::value<int>(diverge)->default_value("0"))
   ("dfac", "cusp divergence exponent for spherical model",
     cxxopts::value<double>(dfac)->default_value("1.0"))
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

				// Get default model bounds unless
				// specificed on the command line
  SphericalModelTable model(filename);

  if (rmin<0.0) rmin = model.get_min_radius();
  if (rmax<0.0) rmax = model.get_max_radius();

				// Generate Sturm-Liouville grid
  auto ortho = std::make_shared<SLGridSph>(filename, Lmax, nmax, numr, rmin, rmax, 
					   true, cmap, rs, 0, 1.0, cachefile, true);
  //                                       ^               ^                  ^
  //                                       |               |                  |
  // Use cache file------------------------+               |                  |
  //                                                       |                  |
  // Cusp extrapolation------------------------------------+                  |
  //                                                                          |
  // Turn on diagnostic output in SL creation---------------------------------+

				// Slaves exit
  if (use_mpi && myid>0) {
    MPI_Finalize();
    exit(0);
  }

				// Do what?
  while (1) {
    bool done=false;
    int iwhich;

    std::cout << "Task:" << std::endl;
    std::cout << "1: Print out density, potential pairs" << std::endl;
    std::cout << "2: Check orthogonality" << std::endl;
    std::cout << "3: Orthogonality matrix" << std::endl;
    std::cout << "4: Quit" << std::endl;
    std::cout << "?? ";
    std::cin >> iwhich;

    switch(iwhich) {
    case 1:
      {
	std::string filename;
	std::cout << "Filename? ";
	std::cin >> filename;
	std::ofstream out (filename.c_str());
	if (!out) {
	  std::cout << "Can't open <" << filename << "> for output" << std::endl;
	  break;
	}

	std::cout << "Number of points? ";
	int num;
	std::cin >> num;

	std::cout << "L, Nmin, Nmax? ";
	int L, Nmin, Nmax;
	std::cin >> L;
	std::cin >> Nmin;
	std::cin >> Nmax;

	Nmin = std::max<int>(Nmin, 0);
	Nmax = std::min<int>(Nmax, nmax);

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
	    << std::setw(15) << " r |";
	for (int n=Nmin; n<Nmax; n++) {
	  std::ostringstream sout1, sout2, sout3;
	  sout1 << " Pot(r, " << n << ") |";
	  sout2 << " Force(r, " << n << ") |";
	  sout3 << " Dens(r, " << n << ") |";
	  out << std::setw(15) << sout1.str()
	      << std::setw(15) << sout2.str()
	      << std::setw(15) << sout3.str();
	}
	out << std::endl;

	out << "# "
	    << std::setw(13) << " [1] |"
	    << std::setw(15) << " [2] |";
	int cntr = 3;
	for (int n=Nmin; n<=Nmax; n++) {
	  std::ostringstream sout1, sout2, sout3;
	  sout1 << "[" << cntr++ << "] |";
	  sout2 << "[" << cntr++ << "] |";
	  sout3 << "[" << cntr++ << "] |";
	  out << std::setw(15) << sout1.str()
	      << std::setw(15) << sout2.str()
	      << std::setw(15) << sout3.str();
	}
	out << std::endl;

	out << "# "
	    << std::setw(13) << std::setfill('-') << "+"
	    << std::setw(15) << std::setfill('-') << "+";
	for (int n=Nmin; n<=Nmax; n++) {
	  out << std::setw(15) << std::setfill('-') << "+"
	      << std::setw(15) << std::setfill('-') << "+"
	      << std::setw(15) << std::setfill('-') << "+";
	}
	out << std::endl << std::setfill(' ');

	// ================
	// END: file header
	// ================

	// =================
	// Begin radial loop
	// =================

	for (int i=0; i<num; i++) {

	  if (use_logr)
	    x = ortho->r_to_xi(exp(lrmin + (lrmax-lrmin)*i/(num-1)));
	  else
	    x = ximin + (ximax-ximin)*i/(num-1);

	  r = ortho->xi_to_r(x);
	  out << setw(15) << x
	      << setw(15) << r;
	  for (int n=Nmin; n<Nmax; n++)
	    out << setw(15) << ortho->get_pot  (x, L, n, 0)
		<< setw(15) << ortho->get_force(x, L, n, 0)
		<< setw(15) << ortho->get_dens (x, L, n, 0);
	  out << std::endl;
	}
      }

      break;

    case 2:
      {
	std::cout << "Number of knots? ";
	int num;
	cin >> num;

	LegeQuad lw(num);

	std::cout << "L, N1, N2? ";
	int L, N1, N2;
	cin >> L;
	cin >> N1;
	cin >> N2;

	double ximin = ortho->r_to_xi(rmin);
	double ximax = ortho->r_to_xi(rmax);

	double x, r, ans=0.0;
	for (int i=0; i<num; i++) {

	  x = ximin + (ximax - ximin)*lw.knot(i);
	  r = ortho->xi_to_r(x);

	  ans += r*r*ortho->get_pot(x, L, N1, 0)*
	    ortho->get_dens(x, L, N2, 0) /
	    ortho->d_xi_to_r(x) * (ximax - ximin)*lw.weight(i);

	}

	std::cout << "<" << N1 << "|" << N2 << "> = " << ans << std::endl;
      }

      break;

    case 3:
      {
	std::cout << "Number of knots? ";
	int num;
	cin >> num;

	LegeQuad lw(num);

	std::cout << "L? ";
	int L;
	std::cin >> L;

	double ximin = ortho->r_to_xi(rmin);
	double ximax = ortho->r_to_xi(rmax);

	Eigen::MatrixXd orthochk(nmax, nmax);

	for (int n1=0; n1<nmax; n1++) {

	  for (int n2=n1; n2<nmax; n2++) {
	    
	    double x, r, ans=0.0;
	    for (int i=0; i<num; i++) {
	  
	      x = ximin + (ximax - ximin)*lw.knot(i);
	      r = ortho->xi_to_r(x);
	      
	      ans += r*r*ortho->get_pot(x, L, n1, 0)*
		ortho->get_dens(x, L, n2, 0) /
		ortho->d_xi_to_r(x) * (ximax - ximin)*lw.weight(i);
	      
	    }
	    
	    orthochk(n1, n2) = orthochk(n2, n1) = ans;
	  }
	}

	std::cout << std::string(60, '-') << std::endl;
	std::cout << orthochk << std::endl;
	std::cout << std::string(60, '-') << std::endl;

      }
      break;

    default:
      done = true;
      break;
    }

    if (done) break;
  }

  if (use_mpi) MPI_Finalize();

  return 0;
}




