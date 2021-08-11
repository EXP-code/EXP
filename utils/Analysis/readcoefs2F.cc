#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "Coefs.H"

int main(int argc, char **argv)
{
  std::string file;
  int nmin, nmax, lmin, lmax;
  bool verbose=false, angle=false, exp_type = true;

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("readcoef",
     "using readcoef output")
    ("nmin",
     po::value<int>(&nmin)->default_value(0), 
     "minimum order for radial coefficients")
    ("nmax",
     po::value<int>(&nmax)->default_value(6), 
     "maximum order for radial coefficients")
    ("lmin",
     po::value<int>(&lmin)->default_value(0), 
     "minimum harmonic order")
    ("lmax",
     po::value<int>(&lmax)->default_value(4), 
     "maximum harmonic order")
    ("file",
     po::value<std::string>(&file)->default_value("coef.dat"),
     "coefficient file")
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  const std::string overview = "Read disk coefficient file and tabulate coefficients for each harmonic subspace in time\n";

  if (vm.count("help")) {
    std::cout << overview << std::endl;
    std::cout << desc     << std::endl;
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
    SphCoefsPtr c = std::make_shared<SphCoefs>();
    if (not c->read(in, exp_type)) break;

    coefs[c->header.tnow] = c;
    std::cout << "Size: " << coefs.size() << std::endl;
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
	    std::cout << std::setw(18) << c.second->coefs[I+s][nn];
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
