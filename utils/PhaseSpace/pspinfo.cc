/*
  Info on each stanza of a phase space dump

  MDWeinberg 03/17/02
*/

#include <stdlib.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <random>
#include <string>
#include <memory>

#include <cxxopts.H>
#include <PSP.H>


int
main(int argc, char *argv[])
{
  bool stats    = false;
  bool timeonly = false;
  bool verbose  = false;
  bool angle    = false;
  std::string new_dir("./");
  int mmin;

  // Option parsing
  //
  cxxopts::Options options(argv[0], "Read a PSP file (either OUT or SPL) and print the metadata\n");

  options.add_options()
    ("h,help", "print this help message")
    ("v,verbose", "print verbose output")
    ("s,stats", "compute and print particle and velocity statistics")
    ("T,time", "print system time info only")
    ("d,dir", "use the provided directory as the data file location",
     cxxopts::value<std::string>(new_dir)->default_value("./"))
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

  if (vm.count("PA")) {
    angle = true;
    mmin = std::max<int>(mmin, 1);
  }

  if (vm.count("time")) {
    timeonly = true;
  }

  if (vm.count("stats")) {
    stats = true;
  }

  // Get trailing arguments
  //
  auto files = vm.unmatched();

  // Sanity check: need at least one file
  //
  if (files.size()==0) {
    std::cout << "You need at least one file" << std::endl
	      << options.help() << std::endl;
    exit(-1);
  }

  // Process the file list
  //
  for (auto file : files) {

    std::shared_ptr<PSP> psp;
    try {
      psp = PSP::getPSP(file, new_dir, verbose);
    }
    catch (const std::exception& e)  {
      std::cout << "pspinfo runtime error: " << e.what() << std::endl;
      exit(-1);
    }
    catch (...) {
      std::cout << "pspinfo unknown error" << std::endl;
      exit(-2);
    }
      
    std::cout << std::string(60, '-')  << std::endl
	      << "Filename: " << file << std::endl;
    psp->PrintSummary(cout, stats, timeonly);
  }

  return 0;
}

