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
  const std::string overview = "View coeffcients from a particular time to test the API";

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v, verbose", "verbose output")
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

  CoefClient coefs(infile);

  // Print the available times
  //
  unsigned cnt = 0;
  for (auto v : coefs.Times()) {
    std::cout << std::setw(16) << v;
    cnt++;
    if (cnt % 5 == 0) std::cout << std::endl;
  }
  if (cnt % 5) std::cout << std::endl;

  // Loop and print coefficients
  //
  while (1) {
    std::cout << "Time? ";

    double time;
    std::cin  >> time;

    auto mat = coefs(time);
    if (mat.rows()) {
      std::cout << std::endl << mat << std::endl << std::endl;
    } else {
      break;
    }

  }

  return(0);
}
