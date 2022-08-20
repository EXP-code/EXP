#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <cxxopts.H>
#include <libvars.H>
#include <CoefsToH5.H>

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

  CoefClient coefs(infile);

  if (false) {
    unsigned cnt = 0;
    for (auto v : coefs.Times()) {
      std::cout << std::setw(16) << v;
      cnt++;
      if (cnt % 5 == 0) std::cout << std::endl;
    }
    if (cnt % 5) std::cout << std::endl;
  }

  coefs.WriteH5Coefs(prefix);

  return(0);
}
