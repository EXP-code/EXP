#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <cxxopts.H>
#include <libvars.H>
#include <Coefficients.H>

int main(int argc, char **argv)
{
  std::string infile, prefix;
  bool verbose = false;

  //
  // Parse Command line
  //
  const std::string overview = "Convert native coefficient file to new HDF5 format";

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v, verbose", "verbose output")
    ("c, cylinder", "assume that coefficients are cylindrical type (for old style)")
    ("s, sphere", "assume that coefficients are spherical type (for old style)")
    ("e, extend", "extend coefficient file rather than write a new file")
    ("i,infile", "input coefficient file",
     cxxopts::value<std::string>(infile)->default_value("coef.dat"))
    ("p,prefix", "prefix for h5 coefficient file",
     cxxopts::value<std::string>(prefix)->default_value("new"))
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

  std::shared_ptr<Coefs::Coefs> coefs;

  // These first two are only needed for converting old-style
  // coefficients (pre-magic-number format).  Otherwise, the last
  // form, the factory member, should deduce the type without
  // specification.
  //
  if (vm.count("cylinder"))	
    coefs = std::make_shared<Coefs::CylCoefs>(infile);
  else if (vm.count("sphere"))
    coefs = std::make_shared<Coefs::SphCoefs>(infile);
  else
    coefs = Coefs::Coefs::factory(infile);

  // Do the writing
  //
  if (vm.count("extend"))
    coefs->ExtendH5Coefs(prefix);
  else
    coefs->WriteH5Coefs(prefix);

  return(0);
}
