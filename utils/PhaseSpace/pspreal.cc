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
  cerr << "Usage: " << prog << ": [-t time -4 -8 -v -h] filename output\n\n";
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -4              convert to float (default)\n";
  cerr << "    -8              convert to double\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool real8 = false, verbose = false;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:48vh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case '4':
      real8 = false;
      break;

    case '8':
      real8 = true;
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

  boost::shared_ptr<ifstream> in;
  boost::shared_ptr<ofstream> out;

  if (optind+1 < argc) {

    in = boost::shared_ptr<ifstream>(new ifstream(argv[optind]));
    if (!in.get()) {
      cerr << "Error opening file <" << argv[optind] << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using input filename: " << argv[optind] << endl;

    
    out = boost::shared_ptr<ofstream>(new ofstream(argv[optind+1]));
    if (!out.get()) {
      cerr << "Error opening file <" << argv[optind+1] << "> for output\n";
      exit(-1);
    }

    if (verbose) cerr << "Using output filename: " << argv[optind+1] << endl;

  } else {
    Usage(prog);
  }


				// Parse the PSP file
				// ------------------
  PSPDump psp(in.get());

				// Now write a summary
				// -------------------
  if (verbose) {

    psp.PrintSummary(in.get(), cerr);
    
    cerr << "\nBest fit dump to <" << time << "> has time <" 
	 << psp.SetTime(time) << ">\n";
  } else 
    psp.SetTime(time);


				// Write the PSP
				// -----------------------------

  in = boost::shared_ptr<ifstream>(new ifstream(argv[optind]));

  psp.writePSP(in.get(), out.get(), !real8);

  return 0;
}
  
