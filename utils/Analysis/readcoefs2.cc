#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <cxxopts.H>
#include "Coefs.H"

int main(int argc, char **argv)
{
  std::string file;
  int nmin, nmax, lmin, lmax;
  bool verbose=false, angle=false, exp_type = true;

  //
  // Parse Command line
  //
  cxxopts::Options options("readcoefs2", "Read disk coefficient file and tabulate coefficients for each harmonic subspace in time");

  options.add_options()
    ("help,h",
     "produce this help message")
    ("PA,p",
     "compute position angle rather than amplitude")
    ("verbose,v",
     "verbose output")
    ("readcoef",
     "using readcoef output")
    ("nmin",
     "minimum order for radial coefficients"
     cxxopts::value<int>(&nmin)->default_value("0"))
    ("nmax",
     "maximum order for radial coefficients",
     cxxopts::value<int>(&nmax)->default_value("6"), 
    ("lmin",
     "minimum harmonic order",
     cxxopts::value<int>(&lmin)->default_value("0"))
    ("lmax",
     "maximum harmonic order",
     cxxopts::value<int>(&lmax)->default_value("4")), 
    ("file",
     "coefficient file",
     cxxopts::value<std::string>(&file)->default_value("coef.dat"))
    ;
  
  auto vm = options.parse(argc, argv);

  const std::string overview = 

  if (vm.count("help")) {
    std::cout << overview << std::endl;
    std::cout << desc     << std::endl;
    return 1;
  }

  if (vm.count("verbose")) verbose = true;

  if (vm.count("PA")) {
    angle = true;
    lmin = std::max<int>(lmin, 1);
  }

  if (vm.count("readcoef")) exp_type = false;

  std::ifstream in(file);
  if (not in) {
    std::cout << "Error opening <" << file << ">" << std::endl;
    return(1);
  }


  std::map<double, SphCoefsPtr> coefs;

  while (in) {
    SphCoefsPtr c = std::make_shared<SphCoefs>();
    if (not c->read(in, exp_type)) break;

    coefs[c->header.tnow] = c;
  }
  
  for (auto c : coefs) {
    unsigned I = 0;
    if (lmin>0) I += lmin*lmin;

    for (int ll=lmin; ll<=std::min<int>(lmax, c.second->header.Lmax); ll++) {
      for (int mm=0; mm<=ll; mm++) {
	std::cout << std::setw(18) << c.first << std::setw(5) << ll << std::setw(5) << mm;
	for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->header.nmax); nn++) {
	  if (mm==0) {
	    if (angle)
	      std::cout << std::setw(18) << 0.0;
	    else
	      std::cout << std::setw(18) << fabs(c.second->coefs(I, nn));
	  } else {
	    if (angle) {
	      double arg = atan2(c.second->coefs(I+1, nn), c.second->coefs(I, nn));
	      std::cout << std::setw(18) << arg;
	    } else {
	      double amp =
		c.second->coefs(I+0, nn) * c.second->coefs(I+0, nn) +
		c.second->coefs(I+1, nn) * c.second->coefs(I+1, nn) ;
	      std::cout << std::setw(18) << sqrt(amp);
	    }
	  }
	}
	std::cout << std::endl;

	if (mm==0) I += 1;
	else       I += 2;

      } // M loop

    } // L loop

  } // T loop
					     

  return(0);
}
