#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>

#include <cxxopts.H>
#include <libvars.H>
#include <CoefContainer.H>

int main(int argc, char **argv)
{
  std::string infile, prefix;
  bool verbose = false;

  //
  // Parse Command line
  //
  const std::string overview = "Generate an ascii table for total harmonic power from the new HDF5 format";

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v, verbose", "verbose output")
    ("i,infile", "coefficient file",
     cxxopts::value<std::string>(infile)->default_value("coef.dat"))
    ("p,prefix", "prefix for the output data file",
     cxxopts::value<std::string>(prefix)->default_value("power"))
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

  auto coefs = Coefs::Coefs::factory(infile);

  auto power = coefs->Power();
  auto times = coefs->Times();

  std::ofstream out(prefix + ".dat");
  if (out) {
    for (int t=0; t<times.size(); t++) {
      out << std::setw(18) << times[t];
      for (int c=0; c<power.cols(); c++)
	out << std::setw(18) << power(t, c);
      out << std::endl;
    }
  }

  return(0);
}
