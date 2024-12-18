#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <random>
#include <cmath>

#include <cxxopts.H>
#include <gaussQ.H>
#include <SLGridMP2.H>
#include "Model1d.H"

int 
main(int argc, char** argv)
{
  double H, zmax, zend;
  int kmax, nmax, knots, numz;
  std::string filename;

  // Parse command line
  //
  cxxopts::Options options(argv[0], "Check the consistency a spherical SL basis");

  options.add_options()
    ("h,help", "Print this help message")
    ("ortho", "Compute orthogonality matrix")
    ("H,height", "Slab scale height",
     cxxopts::value<double>(H)->default_value("1.0"))
    ("K,kmax", "maximum order of in-plane harmonics",
     cxxopts::value<int>(kmax)->default_value("4"))
    ("N,nmax", "maximum number of vertical harmonics",
     cxxopts::value<int>(nmax)->default_value("10"))
    ("n,numz", "size of vertical grid",
     cxxopts::value<int>(numz)->default_value("1000"))
    ("Z,zmax", "maximum extent of vertical grid",
     cxxopts::value<double>(zmax)->default_value("10.0"))
    ("zend",   "potential offset",
     cxxopts::value<double>(zend)->default_value("0.0"))
    ("knots", "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("40"))
    ("p,prefix", "Output filename prefix",
     cxxopts::value<std::string>(filename)->default_value("slabchk_test"))
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

  // Generate Sech2 disk grid
  //
  double dispz = 2.0*M_PI*H*H;

  Sech2  sech2(dispz);
  double h = sech2.get_scale_height();

  SLGridSlab::H = h;
  SLGridSlab::ZEND = zend;
  std::cout << "Check...scale height is: " <<  SLGridSlab::H 
	    << std::endl << std::endl;

  // Particle position generator
  //
  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  auto sample = [&dis, &gen, &h]()
  {
    double m = dis(gen);
    return 0.5*h*log(m/(1.0 - m));
  };

  // Generate Sturm-Liouville grid
  //
  auto ortho = std::make_shared<SLGridSlab>(kmax, nmax, numz, zmax, "isothermal");

  LegeQuad lw(knots);
	  
  double ximin = ortho->z_to_xi(-zmax);
  double ximax = ortho->z_to_xi( zmax);
	  
  std::vector<double> ans1(nmax, 0.0), ans2(nmax, 0.0), ans3(nmax, 0.0);
  
  for (int i=0; i<knots; i++) {
	    
    double x = ximin + (ximax - ximin)*lw.knot(i);
    double z = ortho->xi_to_z(x);

    for (int n=0; n<nmax; n++) {
	      
      ans1[n] += -ortho->get_pot(x, 0, 0, n, 0)*
	ortho->get_dens(x, 0, 0, 1, 0) /
	ortho->d_xi_to_z(x) * (ximax - ximin)*lw.weight(i);
    
      ans2[n] += -ortho->get_pot(x, 0, 0, n, 0)*
	4.0*M_PI*sech2.get_density(z) /
	ortho->d_xi_to_z(x) * (ximax - ximin)*lw.weight(i);
    }
  }
    
  // Monte Carlo version
  //
  int Number = 100000;
  double fac = 4.0*M_PI*2.0*h/Number;

  for (int i=0; i<Number; i++) {
	    
    double z = sample();
    double x = ortho->z_to_xi(z);

    for (int n=0; n<nmax; n++) {
      ans3[n] += -ortho->get_pot(x, 0, 0, n, 0) * fac;
    }
  }
    
  for (int n=0; n<nmax; n++) {
    std::cout << std::setw(6)  << n
	      << std::setw(18) << ans1[n]
	      << std::setw(18) << ans2[n]
	      << std::setw(18) << ans3[n]
	      << std::setw(18) << ans2[n]*ortho->get_dens(0.0, 0, 0, n)/(4.0*M_PI)
	      << std::endl;
  }
  
  std::cout << std::endl;

  int NUM = 20;
  double dz = 3.0*H/(NUM-1);

  for (int i=0; i<NUM; i++) {
	    
    double z = dz*i;
    double x = ortho->z_to_xi(z);
    double s = 0.0, t = 0.0, v = 0.0;

    for (int n=0; n<nmax; n++) {
      s += ans2[n]*ortho->get_dens (x, 0, 0, n, 0);
      t += ans2[n]*ortho->get_pot  (x, 0, 0, n, 0);
      v += ans2[n]*ortho->get_force(x, 0, 0, n, 0);
    }

    std::cout << std::setw(16) << z
	      << std::setw(16) << s/(4.0*M_PI)
	      << std::setw(16) << sech2.get_density(z)
	      << std::setw(16) << t
	      << std::setw(16) << sech2.get_pot(z)
	      << std::setw(16) << v
	      << std::setw(16) << sech2.get_dpot(z)
	      << std::endl;
  }

  std::ofstream out(filename + ".basis");
  if (out) {
    
    int NUM = 200;
    double zmin = -3.0*H, zmax = 3.0*H;
    double dz = (zmax - zmin)/(NUM-1);

    for (int i=0; i<NUM; i++) {
	    
      double z = zmin + dz*i;
      double x = ortho->z_to_xi(z);
      double s = 0.0;

      out << std::setw(16) << z;
      for (int n=0; n<nmax; n++) {
	out << std::setw(16) << ortho->get_pot  (z, 0, 0, n);
	out << std::setw(16) << ortho->get_dens (z, 0, 0, n);
	out << std::setw(16) << ortho->get_force(z, 0, 0, n);
      }
      out << std::endl;
    }
  } else {
    throw std::runtime_error("Error opening filename <" + filename + ".basis>");
  }

  out.close();
  out.open(filename + ".ortho");
  if (out) {
    auto test = ortho->orthoCheck();
    int cnt = 0;
    for ( auto & v : test) {
      out << "==== " << cnt++ << std::endl << v << std::endl;
    }
    
  } else {
    throw std::runtime_error("Error opening filename <" + filename + ".ortho>");
  }


  return 0;
}

