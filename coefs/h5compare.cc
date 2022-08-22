#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <cxxopts.H>
#include <libvars.H>
#include <H5Coefs.H>

int main(int argc, char **argv)
{
  std::string infile, prefix;
  bool verbose = false;

  //
  // Parse Command line
  //
  const std::string overview = "Compare native coefficient file to new HDF5 format";

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v, verbose", "verbose output")
    ("e, extend", "extend coefficient fileb")
    ("i,infile", "native input coefficient file",
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

  CoefClient coefs0(infile);
  CoefClient coefs1(prefix + ".h5");

  // Is data identical?
  //
  if (coefs0.CompareStanzas(coefs1))
    std::cout << "Success" << std::endl;
  else
    std::cout << "FAILURE" << std::endl;


  return(0);
}
