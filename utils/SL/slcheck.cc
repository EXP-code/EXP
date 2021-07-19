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

char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

int main(int argc, char** argv)
{
  bool use_mpi, use_logr;
  double scale, rmin, rmax, rs, dfac;
  int numr, cmap, diverge, Lmax, nmax;
  string filename, cachefile;

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
    ("filename",            po::value<string>(&filename)->default_value("SLGridSph.model"),
     "model file")
    ("cache",               po::value<string>(&cachefile)->default_value(".slgrid_sph_cache"),
     "cache file")
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
  po::variables_map vm;
  
  // Parse command line for control and critical parameters
  //
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error on command line: "
			   << e.what() << std::endl;
    MPI_Finalize();
    return -1;
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
					     true, cmap, rs, true);
  //                                         |               |
  // Use cache file--------------------------+               |
  //                                                         |
  // Turn on diagnostic output in SL creation----------------+

				// Slaves exit
  if (use_mpi && myid>0) {
    MPI_Finalize();
    exit(0);
  }

				// Do what?
  while (1) {
    bool done=false;
    int iwhich;

    cout << "Task:" << endl;
    cout << "1: Print out density, potential pairs" << endl;
    cout << "2: Check orthogonality" << endl;
    cout << "3: Quit" << endl;
    cout << "?? ";
    cin >> iwhich;

    switch(iwhich) {
    case 1:
      {
	std::string filename;
	cout << "Filename? ";
	cin >> filename;
	std::ofstream out (filename.c_str());
	if (!out) {
	  cout << "Can't open <" << filename << "> for output" << endl;
	  break;
	}

	cout << "Number of points? ";
	int num;
	cin >> num;

	cout << "L, Nmin, Nmax? ";
	int L, Nmin, Nmax;
	cin >> L;
	cin >> Nmin;
	cin >> Nmax;

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
	for (int n=Nmin; n<=Nmax; n++) {
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
	  for (int n=Nmin; n<=Nmax; n++)
	    out << setw(15) << ortho->get_pot  (x, L, n, 0)
		<< setw(15) << ortho->get_force(x, L, n, 0)
		<< setw(15) << ortho->get_dens (x, L, n, 0);
	  out << endl;
	}
      }

      break;

    case 2:
      {
	cout << "Number of knots? ";
	int num;
	cin >> num;

	LegeQuad lw(num);

	cout << "L, N1, N2? ";
	int L, N1, N2;
	cin >> L;
	cin >> N1;
	cin >> N2;

	double ximin = ortho->r_to_xi(rmin);
	double ximax = ortho->r_to_xi(rmax);

	double x, r, ans=0.0;
	for (int i=0; i<num; i++) {

	  x = ximin + (ximax - ximin)*lw.knot(i+1);
	  r = ortho->xi_to_r(x);

	  ans += r*r*ortho->get_pot(x, L, N1, 0)*
	    ortho->get_dens(x, L, N2, 0) /
	    ortho->d_xi_to_r(x) * (ximax - ximin)*lw.weight(i+1);

	}

	cout << "<" << N1 << "|" << N2 << "> = " << ans << endl;
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




