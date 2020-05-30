#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/make_unique.hpp>

#include <yaml-cpp/yaml.h>	// YAML support

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
  // Attempt to read coefficient magic number
  //
  const unsigned int cmagic = 0xc0a57a3;
  unsigned int tmagic;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

  if (tmagic == cmagic) {

    // YAML size
    //
    unsigned ssize;
    in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

    // Make and read char buffer
    //
    auto buf = boost::make_unique<char[]>(ssize);
    in.read(buf.get(), ssize);

    YAML::Node node = YAML::Load(buf.get());
      
    // Get parameters
    //
    time = node["time"].as<double>();
    nmax = node["nmax"].as<int>();
    mmax = node["mmax"].as<int>();

  } else {

    // Rewind file
    //
    in.clear();
    in.seekg(0);

    CylCoefHeader header;
    in.read((char *)&header, sizeof(CylCoefHeader));

    time = header.time;
    nmax = header.nmax;
    mmax = header.mmax;
  }

  if (not in) return false;

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
  int nmin, nmax, mmin, mmax;
  bool verbose=false, angle=false;

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
    ("nmin",
     po::value<int>(&nmin)->default_value(0), 
     "minimum order for radial coefficients")
    ("nmax",
     po::value<int>(&nmax)->default_value(6), 
     "maximum order for radial coefficients")
    ("mmin",
     po::value<int>(&mmin)->default_value(0), 
     "minimum azimuthal order")
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

  const std::string overview = "Read disk coefficient file and tabulate coefficients for each harmonic subspace in time\n";

  if (vm.count("help")) {
    std::cout << overview << std::endl;
    std::cout << desc     << std::endl;
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


  std::map<double, CoefPtr> coefs;

  while (in) {
    CoefPtr c = std::make_shared<Coefs>();
    if (not c->read(in)) break;

    coefs[c->time] = c;
  }
  
  for (auto c : coefs) {
    for (int mm=mmin; mm<=std::min<int>(mmax, c.second->mmax); mm++) {
      std::cout << std::setw(18) << c.first << std::setw(5) << mm;
      for (int nn=std::max<int>(nmin, 0); nn<=std::min<int>(nmax, c.second->nmax); nn++) {
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
