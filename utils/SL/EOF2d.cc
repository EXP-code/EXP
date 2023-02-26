#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <PotRZ.H>		// Hankel computation for potential
#include <EmpCyl2D.H>		// 2d empirical basis
#include <cxxopts.H>


int main(int argc, char** argv)
{
  bool logr = false, cmap = false, ortho = false, plane = false;
  int numr, mmax, nmax, knots, M, N, nradial;
  double A, rmin, rmax;
  std::string filename, type, biorth;

  // Parse command line
  //
  cxxopts::Options options(argv[0],
			   "Computes an EOF two-dimensional disk basis from the Clutton-Brock basis for\n"
			   "one of the Kuzmin, finite Mestel or Exponential disk targets.  [The Mestel\n"
			   "disk will work very poorly because the Clutton-Brock basis has infinite\n"
			   "support and looks nothing like the Mestel disk profile.] The new basis, the\n"
			   "orthgonogality matrix and the rotation matrices may be written to files.\n");
  options.add_options()
    ("h,help",     "Print this help message")
    ("logr",       "Plot output grid with logarithmic spacing")
    ("cmap",       "Use mapped coordinates")
    ("ortho",      "Compute EOF orthogonal matrix and write to a file")
    ("grid",       "Print the new basis grid to a file")
    ("trans",      "Print the rotation matrices to a file")
    ("plane",      "Compare the vertical evaluation on the disk plane")
    ("Sk",         "Evaluation the forward transform")
    ("vertical",   "Compute the vertical grid")
    ("basis",      "Use fiducial basis in EmpCyl2D")
    ("rforce",     "Evaluate radial force instead of potential")
    ("zforce",     "Evaluate vertical force instead of potential")
    ("debug",      "Check unitarity in QDHT")
    ("full",       "Use full transform rather than grid evaluation")
    ("totforce",   "Compute the total radial force")
    ("M,harmonic", "Aximuthal harmonic m=0,1,2,3,...",
     cxxopts::value<int>(M)->default_value("0"))
    ("N,norder",   "Default number of knots",
     cxxopts::value<int>(N)->default_value("256"))
    ("n,nradial",   "Radial order for vertical potential output",
     cxxopts::value<int>(nradial)->default_value("0"))
    ("A,length",    "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("mmax",        "maximum number of angular harmonics in the expansion",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("nmax",        "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
    ("numr",        "radial knots for the SL grid",
     cxxopts::value<int>(numr)->default_value("4000"))
    ("r,rmin",      "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
    ("R,rmax",      "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("20.0"))
    ("knots",       "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("200"))
    ("type",        "Target model type (kuzmin, mestel, expon)",
     cxxopts::value<std::string>(type)->default_value("expon"))
    ("biorth",      "Biorthogonal type (cb, bess)",
     cxxopts::value<std::string>(biorth)->default_value("bess"))
    ("o,filename",  "Output filename",
     cxxopts::value<std::string>(filename)->default_value("testeof"))
    ;

  
  //===================
  // Parse options
  //===================

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return 2;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    std::cout << options.help() << std::endl << std::endl;
    return 1;
  }

  // Log spacing?
  //
  if (vm.count("logr")) logr = true;

  // Mapped coordinates?
  //
  if (vm.count("cmap")) cmap = true;

  // Make the class instance
  //
  EmpCyl2D emp(mmax, nmax, knots, numr, rmin, rmax, A, 1.0, cmap, logr,
	       type, biorth);

  if (vm.count("basis")) emp.basisTest(true);
  if (vm.count("debug")) QDHT::debug = true;

  // Sanity check
  //
  M = std::min<int>(M, mmax);

  if (vm.count("grid"))  emp.writeBasis(M, filename + ".grid");
  if (vm.count("trans")) emp.writeTrans(M, filename + ".trans");
  if (vm.count("ortho")) emp.orthoCheck(M, filename + ".ortho");

  emp.checkCoefs();

  if (vm.count("vertical")) {

    // Create the functor
    //
    auto dens = [&emp, M, nradial](double R)
    {
      return emp.get_dens(R, M, nradial);
    };

    // Vertical grid size
    //
    constexpr int num = 40;

    // Output file for grid
    //
    std::ofstream out(filename + ".RZ");

    // Define some representative limits
    //
    double Rmax = 4.0*A;
    double Zmax = 4.0*A;

    // Grid spacing
    //
    double dR = Rmax/(num - 1);
    double dz = Zmax/(num - 1);

    // Get field type
    //
    PotRZ::Field F = PotRZ::Field::potential;
    if (vm.count("rforce")) F = PotRZ::Field::rforce;
    if (vm.count("zforce")) F = PotRZ::Field::zforce;
      
    // Potential instance with radially sensitive convergence parameters
    //
    PotRZ pot(rmax, N, M);

    if (vm.count("full")) {

      out << std::setw(8) << N << std::setw(8) << num << std::endl;

      // Do the grid computation
      //
      for (int j=0; j<num; j++) {
	
	double z = dz*j;
	
	Eigen::VectorXd r(N), p(N);
	
	std::tie(r, p) = pot(z, dens, F);
	
	for (int i=0; i<N; i++) {
	  out << std::setw(16) << r[i]
	      << std::setw(16) << z
	      << std::setw(16) << p[i]
	      << std::endl;
	}
      }
    }
    else if (vm.count("totforce")) {

      out << std::setw(8) << num << std::setw(8) << num << std::endl;

      // Do the grid computation
      //
      for (int j=0; j<num; j++) {
	
	double z = dz*j;
	
	for (int i=0; i<num; i++) {
	  double R = dR*i;
	  double fR = pot(R, z, dens, PotRZ::Field::rforce);
	  double fz = pot(R, z, dens, PotRZ::Field::zforce);
	  out << std::setw(16) << R
	      << std::setw(16) << z
	      << std::setw(16) << sqrt(fR*fR + fz*fz)
	      << std::endl;
	}
      }
    }
    else {

      out << std::setw(8) << num << std::setw(8) << num << std::endl;

      // Do the grid computation
      //
      for (int j=0; j<num; j++) {
	
	double z = dz*j;
	
	for (int i=0; i<num; i++) {
	  out << std::setw(16) << dR*i
	      << std::setw(16) << dz*j
	      << std::setw(16) << pot(dR*i, dz*j, dens, F)
	      << std::endl;
	}
      }
    }
  }

  if (vm.count("plane")) {

    // Output file for grid
    //
    std::ofstream out(filename + ".R0");

    // Potential instance with radially sensitive convergence parameters
    //
    PotRZ pot(rmax, N, M);

    if (vm.count("full")) {

      Eigen::VectorXd r(N), p(N);
      Eigen::MatrixXd outP(nmax, N);

      for (int n=0; n<nmax; n++) {
	// Set the functor using a lambda
	//
	auto dens = [&emp, M, n](double R) { return
	    emp.get_dens(R, M, n);
	};
      
	std::tie(r, p) = pot(0.0, dens);
	outP.row(n) = p;
      }

      // Write the results
      //
      for (int i=0; i<N; i++) {
	out << std::setw(16) << r[i];
	for (int n=0; n<nmax; n++) {
	  out << std::setw(16) <<  outP(n, i)
	      << std::setw(16) << -emp.get_potl(r[i], M, n);
	}
	out << std::endl;
      }

    } else {

      const int ngrid = 40;
      double Rmax = 4.0*A;
      double dR = Rmax/(ngrid-1);

      Eigen::MatrixXd outP(nmax, ngrid);

      for (int n=0; n<nmax; n++) {
	// Set the functor using a lambda
	//
	auto dens = [&emp, M, n](double R) { return
	    emp.get_dens(R, M, n);
	};
      
	for (int j=0; j<ngrid; j++) {
	  outP(n, j) = pot(dR*j, 0.0, dens);
	}
      }

      // Write the results
      //
      for (int i=0; i<ngrid; i++) {
	out << std::setw(16) << dR*i;
	for (int n=0; n<nmax; n++) {
	  out << std::setw(16) <<  outP(n, i)
	      << std::setw(16) << -emp.get_potl(dR*i, M, n);
	}
	out << std::endl;
      }
    }
  }

  if (vm.count("Sk")) {

    // Output file for grid
    //
    std::ofstream out(filename + ".Sk");

    // Potential instance with radially sensitive convergence parameters
    //
    PotRZ pot(rmax, N, M);

    Eigen::MatrixXd outS(N, nmax);
    Eigen::VectorXd k(N), s(N);

    for (int n=0; n<nmax; n++) {
      // Set the functor using a lambda
      //
      auto dens = [&emp, M, n](double R){
	return emp.get_dens(R, M, n); };

      std::tie(k, s) = pot.getKT(dens);
      outS.col(n) = s;
    }


    // Write the results
    //
    for (int i=0; i<k.size(); i++) {
      out << std::setw(16) << k[i];
      for (int n=0; n<nmax; n++) 
	out << std::setw(16) << outS(i, n);
      out << std::endl;
    }
  }

  return 0;
}
