#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <PotRZ.H>		// Hankel computation for potential
#include <EmpCyl2d.H>		// 2d empirical basis
#include <cxxopts.H>


int main(int argc, char** argv)
{
  bool logr = false, cmap = false, ortho = false, plane = false;
  int numr, mmax, nmaxfid, nmax, knots, M, N, nradial, nout;
  double A, rmin, rmax, O, ZN, ZM, rout;
  std::string filename, type, biorth;

  // Parse command line
  //
  cxxopts::Options options(argv[0],
			   "Computes an EOF two-dimensional disk basis from the Clutton-Brock basis for\n"
			   "one of the Kuzmin, finite Mestel, Zang or Exponential disk targets.  [The\n"
			   "Mestel disk will work very poorly because the Clutton-Brock basis has infinite\n"
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
    ("N,nsize",    "Default radial grid size",
     cxxopts::value<int>(N)->default_value("256"))
    ("n,nradial",  "Radial order for vertical potential output",
     cxxopts::value<int>(nradial)->default_value("0"))
    ("A,length",   "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("O,outer",    "outer disk truncation",
     cxxopts::value<double>(O)->default_value("10.0"))
    ("1,TaperN",   "Inner taper exponent",
     cxxopts::value<double>(ZN)->default_value("1.0"))
    ("2,TaperM",   "Outer taper exponent",
     cxxopts::value<double>(ZM)->default_value("4.0"))
    ("mmax",       "maximum number of angular harmonics in the expansion",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("nmaxfid",    "maximum number of radial basis harmonics for EOF construction",
     cxxopts::value<int>(nmaxfid)->default_value("64"))
    ("nmax",       "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
    ("numr",       "radial knots for the SL grid",
     cxxopts::value<int>(numr)->default_value("4000"))
    ("r,rmin",     "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
    ("R,rmax",     "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("20.0"))
    ("knots",      "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("200"))
    ("rout",       "Outer radius for evaluation",
     cxxopts::value<double>(rout)->default_value("10.0"))
    ("nout",       "number of points in the output grid per side",
     cxxopts::value<int>(nout)->default_value("40"))
    ("type",       "Target model type (kuzmin, mestel, zang, expon)",
     cxxopts::value<std::string>(type)->default_value("expon"))
    ("biorth",     "Biorthogonal type (cb, bess)",
     cxxopts::value<std::string>(biorth)->default_value("bess"))
    ("o,filename", "Output filename",
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

  // Parameter vector for the EmpCyl2d models
  //
  EmpCyl2d::ParamVec par {A, O, ZN, ZM};

  // Make the class instance
  //
  EmpCyl2d emp(mmax, nmaxfid, nmax, knots, numr, rmin, rmax, par, A, cmap, logr,
	       type, biorth);

  if (vm.count("basis")) emp.basisTest(true);
  if (vm.count("debug")) QDHT::debug = true;

  // Sanity check
  //
  M = std::min<int>(M, mmax);

  if (vm.count("grid"))  emp.writeBasis(M, filename + ".grid");
  if (vm.count("trans")) emp.writeTrans(M, filename + ".trans");
  if (vm.count("ortho")) emp.orthoCheck(M, filename + ".ortho");

  // Get field type
  //
  PotRZ::Field F = PotRZ::Field::potential;
  if (vm.count("rforce")) F = PotRZ::Field::rforce;
  if (vm.count("zforce")) F = PotRZ::Field::zforce;

  emp.checkCoefs();

  if (vm.count("vertical")) {

    // Create the functor
    //
    auto dens = [&emp, M, nradial](double R)
    {
      return emp.get_dens(R, M, nradial);
    };

    // Output file for grid
    //
    std::ofstream out(filename + ".RZ");

    // Define some representative limits
    //
    double Rmax = rout;
    double Zmax = rout;

    // Grid spacing
    //
    double dR = Rmax/(nout - 1);
    double dz = Zmax/(nout - 1);

    // Potential instance with radially sensitive convergence parameters
    //
    PotRZ pot(rmax, N, M);

    if (vm.count("full")) {

      out << std::setw(8) << N << std::setw(8) << nout << std::endl;

      // Do the grid computation
      //
      for (int j=0; j<nout; j++) {
	
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

      out << std::setw(8) << nout << std::setw(8) << nout << std::endl;

      // Do the grid computation
      //
      for (int j=0; j<nout; j++) {
	
	double z = dz*j;
	
	for (int i=0; i<nout; i++) {
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

      out << std::setw(8) << nout << std::setw(8) << nout << std::endl;

      // Do the grid computation
      //
      for (int j=0; j<nout; j++) {
	
	double z = dz*j;
	
	for (int i=0; i<nout; i++) {
	  out << std::setw(16) << dR*i
	      << std::setw(16) << dz*j;

	  for (int n=0; n<nmax; n++) {
	    auto dens = [&emp, M, n](double R)
	    {
	      return emp.get_dens(R, M, n);
	    };
	    out << std::setw(16) << pot(dR*i, dz*j, dens, F);
	  }
	  out << std::endl;
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
	  out << std::setw(16) << outP(n, i)
	      << std::setw(16) << emp.get_potl(r[i], M, n);
	}
	out << std::endl;
      }

    } else {

      double Rmax = rout;
      double dR = Rmax/(nout-1);

      Eigen::MatrixXd outF(nmax, nout);

      for (int n=0; n<nmax; n++) {
	// Set the functor using a lambda
	//
	auto dens = [&emp, M, n](double R) { return
	    emp.get_dens(R, M, n);
	};
      
	for (int j=0; j<nout; j++) {
	  outF(n, j) = pot(dR*j, 0.0, dens, F);
	}
      }

      // Write the results
      //
      for (int i=0; i<nout; i++) {
	out << std::setw(16) << dR*i;
	for (int n=0; n<nmax; n++) {
	  out << std::setw(16) <<  outF(n, i)
	      << std::setw(16) << -emp.get_potl(dR*i, M, n)
	      << std::setw(16) << -emp.get_dens(dR*i, M, n);
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
