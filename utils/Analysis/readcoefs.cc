#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "../../src/coef.H"

struct Coefs
{
  double time;
  int nmax, mmax;
  std::vector< std::vector<double> > cos_c, sin_c;

  bool read(std::istream& in);
};

typedef std::shared_ptr<Coefs> CoefPtr;

bool Coefs::read(std::istream& in)
{
  CoefHeader2 header;
  in.read((char *)&header, sizeof(CoefHeader2));
  if (not in) return false;

  time = header.time;
  nmax = header.nmax;
  mmax = header.mmax;

  cos_c.resize(mmax+1);
  sin_c.resize(mmax+1);
  
  for (int mm=0; mm<=mmax; mm++) {

    cos_c[mm].resize(nmax);
    in.read((char *)&cos_c[mm][0], sizeof(double)*nmax);

    if (mm) {
      sin_c[mm].resize(nmax);
      in.read((char *)&sin_c[mm][0], sizeof(double)*nmax);
    }
  }

  return true;
}


int main(int argc, char **argv)
{
  std::string file;
  int nmax, mmax;
  bool verbose;

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("nmax",
     po::value<int>(&nmax)->default_value(6), 
     "maximum order for radial coefficients")
    ("mmax",
     po::value<int>(&mmax)->default_value(4), 
     "maximum azimuthal order")
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

  if (vm.count("verbose")) verbose = true;

  std::ifstream in(file);
  if (not in) {
    std::cout << "Error opening <" << file << ">" << std::endl;
    return(1);
  }


  std::map<double, CoefPtr> coefs;

  while (in) {
    CoefPtr c = std::make_shared<Coefs>();
    if (not c->read(in)) break;

    coefs[c->time] = c;
  }
  
  for (auto c : coefs) {
    for (int mm=0; mm<=std::min<int>(mmax, c.second->mmax); mm++) {
      std::cout << std::setw(12) << c.first << std::setw(5) << mm;
      for (int nn=0; nn<=std::min<int>(nmax, c.second->nmax); nn++) {
	if (mm==0)
	  std::cout << std::setw(12) << fabs(c.second->cos_c[mm][nn]);
	else {
	  double amp =
	    c.second->cos_c[mm][nn] * c.second->cos_c[mm][nn] +
	    c.second->sin_c[mm][nn] * c.second->sin_c[mm][nn] ;
	  std::cout << std::setw(12) << sqrt(amp);
	}
      }
      std::cout << std::endl;
    }
  }
					     

  return(0);
}
