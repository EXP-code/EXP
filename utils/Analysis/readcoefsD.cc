#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <cxxopts.H>
#include "Coefs.H"

int main(int argc, char **argv)
{
  std::string file;
  int nmin, nmax, mmin, mmax;
  bool verbose=false, angle=false;

  //
  // Parse Command line
  //
  cxxopts::Options options("readcoefs", "Read disk coefficient file and tabulate coefficients for each harmonic subspace in time");

  options.add_options()
    ("h,help",
     "produce this help message")
    ("p,PA",
     "compute position angle rather than amplitude")
    ("v,verbose",
     "verbose output")
    ("nmin",
     "minimum order for radial coefficients",
     cxxopts::value<int>(nmin)->default_value("0"))
    ("nmax",
     "maximum order for radial coefficients",
     cxxopts::value<int>(nmax)->default_value("6"))
    ("mmin",
     "minimum azimuthal order",
     cxxopts::value<int>(mmin)->default_value("0"))
    ("mmax",
     "maximum azimuthal order",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("file",
     "coefficient file",
     cxxopts::value<std::string>(file)->default_value("coef.dat"))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  if (vm.count("verbose")) verbose = true;

  if (vm.count("PA")) {
    angle = true;
    mmin = std::max<int>(mmin, 1);
  }

  std::ifstream in(file);
  if (not in) {
    std::cout << "Error opening <" << file << ">" << std::endl;
    return(1);
  }


  std::map<double, CylCoefsPtr> coefs;

  while (in) {
    CylCoefsPtr c = std::make_shared<CylCoefs>();
    if (not c->read(in, verbose)) break;

    coefs[c->time] = c;
  }
  
  for (auto c : coefs) {
    for (int mm=mmin; mm<=std::min<int>(mmax, c.second->mmax); mm++) {
      std::cout << std::setw(18) << c.first << std::setw(5) << mm;
      for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->nmax); nn++) {
	if (mm==0) {
	  if (angle)
	    std::cout << std::setw(18) << 0.0;
	  else
	    std::cout << std::setw(18) << fabs(c.second->cos_c[mm][nn]);
	} else {
	  if (angle) {
	    double arg = atan2(c.second->sin_c[mm][nn], c.second->cos_c[mm][nn]);
	    std::cout << std::setw(18) << arg;
	  } else {
	    double amp =
	      c.second->cos_c[mm][nn] * c.second->cos_c[mm][nn] +
	      c.second->sin_c[mm][nn] * c.second->sin_c[mm][nn] ;
	    std::cout << std::setw(18) << sqrt(amp);
	  }
	}
      }
      std::cout << std::endl;
    }
  }
					     

  return(0);
}
