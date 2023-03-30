/*
  Test Disk2d by computing coefficients and evaluating fields for an
  exponential disk target
*/
                                // C++/STL headers
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <random>
#include <vector>
#include <cmath>

#include <fenv.h>

#include <config_exp.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

// EXP classes
//
#include <numerical.H>
#include <gaussQ.H>
#include <interp.H>
#include <Disk2d.H>
#include <libvars.H>		// Library globals
#include <cxxopts.H>		// Command-line parsing
#include <EXPini.H>		// Ini-style config

#include <norminv.H>

using namespace __EXP__;	// Reference to n-body globals

int 
main(int ac, char **av)
{
  //====================
  // Inialize MPI stuff
  //====================
  
  local_init_mpi(ac, av);

  //====================
  // Begin opt parsing
  //====================

  double       RCYLMIN, RCYLMAX, DMFAC, ACYL;
  bool         LOGR, CHEBY, SELECT, DUMPCOEF;
  int          NUMR, CMAPR, CMAPZ, knots, nlim, M;
  int          NMAXH, NMAXD, NMAXFID, LMAX, MMAX, NUMX, NUMY, NOUT;
  double       scale_length, disk_mass, mamp;
  std::string  cachefile, config, prefix;
  
  const std::string mesg("Generates a 2d cylindrical basis and computes the fields for an exponential disk target for checking convergence\n");

  cxxopts::Options options(av[0], mesg);

  options.add_options()
    ("h,help", "Print this help message")
    ("T,template", "Write template options file with current and all default values")
    ("c,config", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("cachefile", "The cache file for the cylindrical basis",
     cxxopts::value<string>(cachefile)->default_value(".eof_2d_cache"))
    ("prefix", "Output file prefix",
     cxxopts::value<string>(prefix)->default_value("test2d"))
    ("NOUT", "Number of test output points",
     cxxopts::value<int>(NOUT)->default_value("200"))
    ("NUMR", "Number of grid points in BiorthCyl grid",
     cxxopts::value<int>(NUMR)->default_value("2000"))
    ("NUMX", "Number of knots in radial dimension of meridional grid",
     cxxopts::value<int>(NUMX)->default_value("256"))
    ("NUMY", "Number of knots in vertical dimension of meridional grid",
     cxxopts::value<int>(NUMY)->default_value("128"))
    ("MMAX", "Harmonic order for disk expansion",
     cxxopts::value<int>(MMAX)->default_value("12"))
    ("NMAX", "Number of radial terms in cylindrical expansion",
     cxxopts::value<int>(NMAXD)->default_value("12"))
    ("NMAXFID", "Number of radial terms in Bessel EOF expansion",
     cxxopts::value<int>(NMAXFID)->default_value("64"))
    ("LOGR", "Logarithmic scaling for model grid",
     cxxopts::value<bool>(LOGR)->default_value("true"))
    ("disk_mass", "Mass of the stellar disk",
     cxxopts::value<double>(disk_mass)->default_value("0.05"))
    ("scale_length", "Scale length for the realized disk",
     cxxopts::value<double>(scale_length)->default_value("0.01"))
    ("RCYLMIN", "Minimum disk radius",
     cxxopts::value<double>(RCYLMIN)->default_value("0.001"))
    ("RCYLMAX", "Maximum disk radius, in units of ASCALE",
     cxxopts::value<double>(RCYLMAX)->default_value("20.0"))
    ("ACYL", "Radial scale length for disk basis construction",
     cxxopts::value<double>(ACYL)->default_value("0.6"))
    ("CMAPR", "Radial coordinate mapping type for cylindrical grid  (0=none, 1=rational fct)",
     cxxopts::value<int>(CMAPR)->default_value("1"))
    ("CMAPZ", "Vertical coordinate mapping type for cylindrical grid  (0=none, 1=rational fct)",
     cxxopts::value<int>(CMAPZ)->default_value("1"))
    ("knots", "Knots for scalar product",
     cxxopts::value<int>(knots)->default_value("120"))
    ("nlim", "Limit for the number of radial terms",
     cxxopts::value<int>(nlim)->default_value("1000000"))
    ("report", "Print out progress in BiorthCyl table evaluation")
    ("coefs", "Print out coefficient array for m=0")
    ("M,MP", "Perturbation harmonic",
     cxxopts::value<int>(M)->default_value("2"))
    ("pert", "Harmonic amplitude",
     cxxopts::value<double>(mamp)->default_value("0.0"))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Write YAML template config file and exit
  //
  if (vm.count("template")) {
    // Write template file
    //
    if (myid==0) SaveConfig(vm, options, "template.yaml");

    MPI_Finalize();
    return 0;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << options.help() << std::endl << std::endl
		<< "Examples: " << std::endl
		<< "\t" << "Use parameters read from a YAML config file"  << std::endl
		<< "\t" << av[0] << " --config=gendisk.config"  << std::endl << std::endl
		<< "\t" << "Generate a template YAML config file from current defaults"  << std::endl
		<< "\t" << av[0] << " --template=template.yaml" << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the YAML config file
  //
  if (vm.count("config")) {
    try {
      vm = LoadConfig(options, config);
    } catch (cxxopts::OptionException& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }
  
  // Create YAML db
  //
  YAML::Emitter yml;
  yml << YAML::BeginMap;
  yml << YAML::Key << "acyltbl"   << YAML::Value << ACYL;
  yml << YAML::Key << "rcylmin"   << YAML::Value << RCYLMIN;
  yml << YAML::Key << "rcylmax"   << YAML::Value << RCYLMAX;
  yml << YAML::Key << "scale"     << YAML::Value << scale_length;
  yml << YAML::Key << "numx"      << YAML::Value << NUMX;
  yml << YAML::Key << "numy"      << YAML::Value << NUMY;
  yml << YAML::Key << "numr"      << YAML::Value << NUMR;
  yml << YAML::Key << "cmapR"     << YAML::Value << CMAPR;
  yml << YAML::Key << "cmapZ"     << YAML::Value << CMAPZ;
  yml << YAML::Key << "Mmax"      << YAML::Value << MMAX;
  yml << YAML::Key << "nmax"      << YAML::Value << NMAXD;
  yml << YAML::Key << "nmaxfid"   << YAML::Value << NMAXFID;
  yml << YAML::Key << "logr"      << YAML::Value << LOGR;
  yml << YAML::Key << "cachename" << YAML::Value << cachefile;
  if (vm.count("report"))
    yml << YAML::Key << "verbose" << YAML::Value << true;
  yml << YAML::EndMap;
  
  // Create expansion in parallel if needed . . .
  //
  Disk2d disk(YAML::Load(yml.c_str()));

  // Only root does comparison and output
  //
  if (myid==0) {

    auto dens = [&scale_length, &disk_mass, &mamp, &M](int m, double r)
    {
      double dens0 = disk_mass/(2.0*M_PI*scale_length*scale_length)*exp(-r/scale_length);
      if (m==0) return dens0;
      if (m==M) return dens0*mamp*
		  exp(-(r - scale_length)*(r - scale_length)/(2.0*scale_length*scale_length*0.09));
      else      return 0.0;
    };

    disk.inner_product(dens, knots, nlim);

    std::ofstream out(prefix + ".data");
    if (out) {
      out << "#" << std::setw(15) << "radius |"
	  << std::setw(16) << "dens(m=0) |"
	  << std::setw(16) << "dens(m=2) |"
	  << std::setw(16) << "rho (m=0) |"
	  << std::setw(16) << "rho (m>0) |"
	  << std::setw(16) << "pot (m=0) |"
	  << std::setw(16) << "pot (m>0) |"
	  << std::setw(16) << "Fr |" << std::endl
	  << "#" << std::setw(15) << "[1] |"
	  << std::setw(16) << "[2] |"
	  << std::setw(16) << "[3] |"
	  << std::setw(16) << "[4] |"
	  << std::setw(16) << "[5] |"
	  << std::setw(16) << "[6] |"
	  << std::setw(16) << "[7] |"
	  << std::setw(16) << "[8] |" << std::endl;

      double rmin = log(ACYL*1.0e-4), rmax = log(RCYLMAX);
      double dr = (rmax - rmin)/NOUT;
      double d0, d1, p0, p1, Fr, Fz, Fp;
      for (int i=0; i<=NOUT; i++) {
	double r = exp(rmin + dr*i);
	std::tie(d0, d1, p0, p1, Fr, Fz, Fp) = disk.accumulated_eval(r, 0.0, 0.0);
	
	out << std::setw(16) << r
	    << std::setw(16) << dens(0, r)
	    << std::setw(16) << dens(2, r)
	    << std::setw(16) << d0
	    << std::setw(16) << d1
	    << std::setw(16) << p0
	    << std::setw(16) << p1
	    << std::setw(16) << Fr
	    << std::endl;
      }

      if (vm.count("coefs")) {
	std::cout << "Coefficients m=0: " << disk.get_coefs().row(0)
		  << std::endl;
	std::cout << "Coefficients m=" << M << ": "
		  << disk.get_coefs().row(1+2*(M-1)) << std::endl;
      }
    } else {
      std::cout << av[0] << ": could not open <" << prefix + ".data"
		<< "> for output" << std::endl;
    }
  }

  // Done
  //
  MPI_Finalize();

  return 0;
}
