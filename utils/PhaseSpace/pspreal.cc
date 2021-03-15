/*
  Convert from real quanities to float or double

  MDWeinberg 01/07/13
*/

using namespace std;

#include <cstdlib>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>

#include <StringTok.H>
#include <header.H>
#include <PSP.H>

				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << "Synposis: convert PSP file real fields to " << sizeof(float)
       << "-byte (single) or " << sizeof(double) 
       << "-byte (double) float types" << std::endl << std::endl;
  cerr << "Usage: " << prog << ": [-4 -8 -v -h] filename output\n\n";
  cerr << "    -4              convert to float (default)\n";
  cerr << "    -8              convert to double\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -S              use SPL format\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  bool real8 = false, verbose = false, SPL = false;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:48vh");

    if (c == -1) break;

    switch (c) {

    case '4':
      real8 = false;
      break;

    case '8':
      real8 = true;
      break;

    case 'S':
      SPL = true;
      break;

    case 'v':
      verbose = true;
      break;

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

  std::string file = argv[optind];
  std::string ofil = argv[optind+1];
  std::ofstream out;

  if (optind+1 < argc) {

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
    
  } else {
    Usage(prog);
  }


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
  
