/*
  Convert from real quanities to float or double

  MDWeinberg 01/07/13
*/

#include <cstdlib>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <vector>
#include <string>
#include <list>

#include <StringTok.H>
#include <header.H>
#include <PSP.H>
#include <cxxopts.H>

//-------------
// Help message
//-------------



int
main(int argc, char **argv)
{
  char *prog = argv[0];
  bool real8 = false, verbose = false, SPL = false;
  std::string file, ofil;

  // Parse command line
  //
  cxxopts::Options options("pspreal", "Convert the PSP to and from float and double");
  options.add_options()
    ("h,help", "this help message")
    ("v,verbose", "verbose output")
    ("4,tofloat", "convert to default (default)")
    ("8,todouble", "convert to double")
    ("i,input", "input file", cxxopts::value<std::string>(file))
    ("o,output", "output file", cxxopts::value<std::string>(ofil))
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
    exit(-1);
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  if (vm.count("float")) {
    real8 = false;
  }

  if (vm.count("double")) {
    real8 = true;
  }

  if (vm.count("SPL")) {
    SPL = true;
  }

  if (vm.count("input")==0) {
    std::cout << "You must specify the input file" << std::endl;
    exit(-1);
  }

  if (vm.count("output")==0) {
    std::cout << "You must specify the output file" << std::endl;
    exit(-1);
  }

  std::ofstream out;

  std::ifstream in(file);
  if (!in) {
    std::cerr << "Error opening file <" << file << "> for input" << std::endl;
    exit(-1);
  }

  if (verbose) std::cerr << "Using input filename: " << file << std::endl;

    
  out.open(ofil);
  if (!out) {
    std::cerr << "Error opening file <" << ofil << "> for output" << std::endl;
    exit(-1);
  }

  if (verbose) cerr << "Using output filename: " << ofil << endl;
    

				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (SPL) psp = std::make_shared<PSPspl>(file);
    else     psp = std::make_shared<PSPout>(file);


				// Now write a summary
				// -------------------
  if (verbose) {

    psp->PrintSummary(cerr);
    
    std::cerr << "\nBest fit dump to <" << time << "> has time <" 
	      << psp->CurrentTime() << ">" << std::endl;
  }


				// Write the PSP
				// -----------------------------

  psp->writePSP(out, !real8);

  return 0;
}
  
