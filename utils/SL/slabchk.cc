#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <random>
#include <cmath>

#include "cxxopts.H"
#include "gaussQ.H"
#include "SLGridMP2.H"
#include "Model1d.H"
#include "numerical.H"

int 
main(int argc, char** argv)
{
  double H, zmax, zend, sigma, eps;
  int kx, nmax, knots, numz, nfid;
  std::string model, filename;

  // Parse command line
  //
  cxxopts::Options options(argv[0], "Check the consistency a slab SL basis");

  options.add_options()
    ("h,help", "Print this help message")
    ("o,ortho", "Compute orthogonality matrix")
    ("H,height", "Slab scale height",
     cxxopts::value<double>(H)->default_value("1.0"))
    ("K,kx", "order of in-plane harmonics",
     cxxopts::value<int>(kx)->default_value("0"))
    ("N,nmax", "maximum number of vertical harmonics",
     cxxopts::value<int>(nmax)->default_value("18"))
    ("n,numz", "size of vertical grid",
     cxxopts::value<int>(numz)->default_value("1000"))
    ("Z,zmax", "maximum extent of vertical grid",
     cxxopts::value<double>(zmax)->default_value("10.0"))
    ("s,sigma","Velocity dispersion for model",
     cxxopts::value<double>(sigma)->default_value("1.0"))
    ("d,zend",  "potential offset",
     cxxopts::value<double>(zend)->default_value("0.0"))
    ("e,eps",   "edge offset for basis and sampling",
     cxxopts::value<double>(eps)->default_value("1.0e-3"))
    ("k,knots", "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("40"))
    ("j,ncheck", "Order of coefficient to check for orthogonality",
     cxxopts::value<int>(nfid)->default_value("2"))
    ("m,model", "Model type (Sech2mu, Uniform)",
     cxxopts::value<std::string>(model)->default_value("Sech2mu"))
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

  // Define model type and instantiate model
  //
  OneDModel::ModelType mtype = OneDModel::parse_model_type(model);
  std::shared_ptr<OneDModel> mdl;
  std::string mname = "isothermal";
  
  if (mtype == OneDModel::ModelType::uniform) {
    
    std::cout << "Harmonic potential: setting scale height to "
	      << H << " and Zmax to " << H*(1.0 - eps) << std::endl;
    mdl = std::make_shared<Uniform>(H, sigma, sigma);
    zmax  = H*(1.0 - eps);
    mname = "uniform";
  }
  else if (mtype == OneDModel::ModelType::cosine) {
    
    std::cout << "Cosine bell potential: setting scale height to "
	      << H << " and Zmax to " << H*0.999 << std::endl;
    mdl = std::make_shared<Cosine>(H, sigma, sigma);
    zmax  = H*0.999;
    mname = "cosine";
  } else if (mtype == OneDModel::ModelType::sech2mu) {
    
    std::cout << "Sech2 potential: using scale height " << H << std::endl;

    auto tmp = std::make_shared<Sech2mu>(sigma*sigma, H);
    tmp->set_hmax(10.0*H);
    mdl = tmp;
    
    SLGridSlab::H = H;
    SLGridSlab::ZEND = zend;
    mname = "sech2mu";
  } else {
    std::cerr << "model must be 'Sech2mu' or 'Harmonic'\n";
    return 1;
  }

  // Particle position generator
  //
  std::random_device rd; 
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  auto sample = [&dis, &gen, &H, &eps, &mtype]()
  {
    if (mtype == OneDModel::ModelType::uniform) {
      return H*(2.0*dis(gen) - 1.0);
    }
    else if (mtype == OneDModel::ModelType::cosine) {
      // The cumulative mass profile for the cosine bell model is:
      //
      auto M = [&H](double z) -> double {
	return 0.5/H*(z + H + H/M_PI*sin(M_PI*z/H));
      };

      // Generate a random mass and invert M(z) using Brent's method to get z
      //
      double del  = 1.0 - eps;
      double Mbot = M(-H*del), Mtop = M(H*del);
      double m = Mbot + (Mtop - Mbot)*dis(gen);

      auto R = [&m, &H](double z) -> double {
	return 0.5/H*(z + H + H/M_PI*sin(M_PI*z/H)) - m;
      };

      return zbrent(R, -H*del, H*del, 1e-10);
    }
    else if (mtype == OneDModel::ModelType::sech2mu) {
      double m = dis(gen);
      return 0.5*H*log(m/(1.0 - m));
    }
    else {
      throw std::runtime_error("Invalid model type for sampling");
    }
  };

  // Generate Sturm-Liouville grid
  //
  auto ortho = std::make_shared<SLGridSlab>(kx+1, nmax, numz, zmax,
					    ".slgrid_slab_cache",
					    "linear", mname, false);

  LegeQuad lw(knots);
	  
  double ximin = ortho->z_to_xi(-zmax);
  double ximax = ortho->z_to_xi( zmax);
	  
  std::vector<double> ans1(nmax, 0.0), ans2(nmax, 0.0), ans3(nmax, 0.0);
  
  for (int i=0; i<knots; i++) {
	    
    double x = ximin + (ximax - ximin)*lw.knot(i);
    double z = ortho->xi_to_z(x);

    for (int n=0; n<nmax; n++) {
	      
      ans1[n] += -ortho->get_pot(x, kx, 0, n, 0)*
	ortho->get_dens(x, 0, 0, nfid, 0) /
	ortho->d_xi_to_z(x) * (ximax - ximin)*lw.weight(i);
    
      ans2[n] += -ortho->get_pot(x, kx, 0, n, 0)*
	4.0*M_PI*mdl->get_density(z) /
	ortho->d_xi_to_z(x) * (ximax - ximin)*lw.weight(i);
    }
  }
    
  // Monte Carlo version
  //
  int Number = 100000;
  double fac = 4.0*M_PI/Number;

  for (int i=0; i<Number; i++) {
	    
    double z = sample();
    double x = ortho->z_to_xi(z);

    for (int n=0; n<nmax; n++) {
      ans3[n] += -ortho->get_pot(x, kx, 0, n, 0) * fac;
    }
  }
    
  std::cout << "---- Coefficients\n";
  std::cout << std::setw(6)  << "n"
	    << std::setw(18) << "Dens(1)"
	    << std::setw(18) << "Dens(model)"
	    << std::setw(18) << "Dens(MC)"
	    << std::setw(18) << "Model dens at z=0"
	    << std::endl;

  for (int n=0; n<nmax; n++) {
    std::cout << std::setw(6)  << n
	      << std::setw(18) << ans1[n]
	      << std::setw(18) << ans2[n]
	      << std::setw(18) << ans3[n]
	      << std::setw(18) << ans2[n]*ortho->get_dens(0.0, kx, 0, n)/(4.0*M_PI)
	      << std::endl;
  }
  
  std::cout << std::endl;

  int NUM = 20;
  double dz = zmax/(NUM-1);

  std::cout << "---- Comparing SL solution to model . . .\n";
  std::cout << std::setw(16) << "z"
	    << std::setw(16) << "rho (SL)"
	    << std::setw(16) << "rho (mod)"
	    << std::setw(16) << "phi (SL)"
	    << std::setw(16) << "phi (mod)"
	    << std::setw(16) << "dphi (SL)"
	    << std::setw(16) << "dphi (mod)"
	    << std::endl;

  for (int i=0; i<NUM; i++) {
	    
    double z = dz*i;
    double x = ortho->z_to_xi(z);
    double s = 0.0, t = 0.0, v = 0.0;

    for (int n=0; n<nmax; n++) {
      s += ans2[n]*ortho->get_dens (x, kx, 0, n, 0);
      t += ans2[n]*ortho->get_pot  (x, kx, 0, n, 0);
      v += ans2[n]*ortho->get_force(x, kx, 0, n, 0);
    }

    std::cout << std::setw(16) << z
	      << std::setw(16) << s/(4.0*M_PI)
	      << std::setw(16) << mdl->get_density(z)
	      << std::setw(16) << t
	      << std::setw(16) << mdl->get_pot(z)
	      << std::setw(16) << v
	      << std::setw(16) << mdl->get_dpot(z)
	      << std::endl;
  }

  std::ofstream out(filename + ".basis");
  if (out) {
    
    int NUM = 200;
    double Zmin = -zmax, Zmax = zmax;
    double dz = (Zmax - Zmin)/(NUM-1);

    for (int i=0; i<NUM; i++) {
	    
      double z = Zmin + dz*i;
      double x = ortho->z_to_xi(z);
      double s = 0.0;

      out << std::setw(16) << z;
      for (int n=0; n<nmax; n++) {
	out << std::setw(16) << ortho->get_pot  (z, kx, 0, n);
	out << std::setw(16) << ortho->get_dens (z, kx, 0, n);
	out << std::setw(16) << ortho->get_force(z, kx, 0, n);
      }
      out << std::endl;
    }
  } else {
    throw std::runtime_error("Error opening filename <" + filename + ".basis>");
  }

  out.close();

  if (vm.count("ortho")) {
    std::cout << "---- Computing orthogonality matrix\n";

    out.open(filename + ".ortho");
    if (out) {
      auto test = ortho->orthoCheck();
      int cnt = 0;
      for ( auto & v : test) {
	out << "==== " << cnt++ << std::endl << v << std::endl;
      }
      
      std::cout << "---- Orthogonality matrix written to <" << filename << ".ortho>\n";

    } else {
      throw std::runtime_error("Error opening filename <" + filename + ".ortho>");
    }
  }
  
  return 0;
}

