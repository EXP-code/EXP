#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include "cxxopts.H"
#include "libvars.H"
#include "Coefs.H"

int main(int argc, char **argv)
{
  std::string file;
  int nmin, nmax, lmin, lmax;
  bool verbose=false, angle=false, exp_type = true;

  //
  // Parse Command line
  //
  const std::string overview = "Read disk coefficient file and tabulate coefficients for each harmonic subspace in time";

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v, verbose", "verbose output")
    ("readcoef", "using readcoef output")
    ("nmin", "minimum order for radial coefficients",
     cxxopts::value<int>(nmin)->default_value("0"))
    ("nmax", "maximum order for radial coefficientns",
     cxxopts::value<int>(nmax)->default_value("6"))
    ("lmin", "minimum harmonic order",
     cxxopts::value<int>(lmin)->default_value("0"))
    ("lmax", "maximum harmonic order",
     cxxopts::value<int>(lmax)->default_value("4"))
    ("file", "coefficient file",
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

  if (vm.count("readcoef")) exp_type = false;

  std::ifstream in(file);
  if (not in) {
    std::cout << "Error opening <" << file << ">" << std::endl;
    return(1);
  }


  std::map<double, SphCoefsPtr> coefs;

  while (in) {
    try {
      SphCoefsPtr c = std::make_shared<SphCoefs>();
      if (not c->read(in, exp_type)) break;

      coefs[c->header.tnow] = c;
    }
    catch(std::runtime_error& error) {
      std::cout << "Error (runtime): " << error.what() << std::endl;
      break;
    }
    catch(std::logic_error& error) {
      std::cout << "Error (logic): " << error.what() << std::endl;
      break;
    }
      
  }
  
  for (auto c : coefs) {
    unsigned I = 0;
    if (lmin>0) I += lmin*lmin;

    for (int ll=lmin; ll<=std::min<int>(lmax, c.second->header.Lmax); ll++) {
      for (int mm=0; mm<=ll; mm++) {
	int S = mm==0 ? 1 : 2;
	for (int s=0; s<S; s++) {
	  std::cout << std::setw(18) << c.first << std::setw(5) << ll << std::setw(5) << mm << std::setw(5) << s;
	  for (int nn=std::max<int>(nmin, 0); nn<std::min<int>(nmax, c.second->header.nmax); nn++) 
	    std::cout << std::setw(18) << c.second->coefs(I+s, nn);
	  std::cout << std::endl;
	  if (mm==0) break;
	}
	// Cosine/sine

	if (mm==0) I += 1;
	else       I += 2;

      } // M loop

    } // L loop

  } // T loop
					     

  return(0);
}
