#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "../../src/coef.H"

struct SphCoefs
{
  SphCoefHeader header;
  std::vector< std::vector<double> > coefs;

  bool read(std::istream& in, bool exp_type);
};

typedef std::shared_ptr<SphCoefs> SphCoefsPtr;

bool SphCoefs::read(std::istream& in, bool exp_type)
{
  in.read((char *)&header, sizeof(SphCoefHeader));
  if (not in) return false;

  coefs.resize((header.Lmax+1)*(header.Lmax+1));
  for (auto & v : coefs) v.resize(header.nmax);

  for (int ir=0; ir<header.nmax; ir++) {
    for (int l=0; l<(header.Lmax+1)*(header.Lmax+1); l++)
      in.read((char *)&coefs[l][ir], sizeof(double));
  }

  if (exp_type) {
    int k = 0;
    for (int l=0; l<=header.Lmax; l++) {
      for (int m=0; m<=l; m++) {
	double fac = sqrt( (0.5*l+0.25)/M_PI * 
			   exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );

	if (m != 0) fac *= M_SQRT2;

	// Cosine terms
	for (int ir=0; ir<header.nmax; ir++) coefs[k][ir] *= fac;
	k++;

	// Sine terms
	if (m != 0) {
	  for (int ir=0; ir<header.nmax; ir++) coefs[k][ir] *= fac;
	  k++;
	}
      }
    }
  }

  return true;
}


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
    ("PA,p",
     "compute position angle rather than amplitude")
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
	for (int nn=std::max<int>(nmin, 0); nn<=std::min<int>(nmax, c.second->header.nmax); nn++) {
	  if (mm==0) {
	    if (angle)
	      std::cout << std::setw(18) << 0.0;
	    else
	      std::cout << std::setw(18) << fabs(c.second->coefs[I][nn]);
	  } else {
	    if (angle) {
	      double arg = atan2(c.second->coefs[I+1][nn], c.second->coefs[I][nn]);
	      std::cout << std::setw(18) << arg;
	    } else {
	      double amp =
		c.second->coefs[I+0][nn] * c.second->coefs[I+0][nn] +
		c.second->coefs[I+1][nn] * c.second->coefs[I+1][nn] ;
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
