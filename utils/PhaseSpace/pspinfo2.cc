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
#include <string>
#include <memory>

#include <PSP2.H>

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
  cerr << prog << ": [-v -t] filename\n";
  cerr << "        -v    verbose output\n";
  cerr << "        -S    assume split PSP files\n";
  cerr << "        -s    particle and velocity statistics\n";
  cerr << "        -T    print time info only\n";
  exit(-1);
}

int
main(int argc, char *argv[])
{
  bool spl = false;
  bool stats = false;
  bool timeonly = false;
  bool verbose = false;
  int c;
  
  while (1) {
    c = getopt(argc, argv, "vsSTh");
    if (c == -1) break;

    switch (c) {
      
    case 'v':
      verbose = true;
      break;

    case 'S':
      spl = true;
      break;

    case 's':
      stats = true;
      break;

    case 'T':
      timeonly = true;
      break;

    case 'h':
    case '?':
    default:
      Usage(argv[0]); 
    
    }
  }

  if (optind >= argc) Usage(argv[0]);

  std::string file(argv[optind]);

  cerr << "Filename: " << file << endl;

  if (not spl and file.find("SPL")!=std::string::npos) spl = true;

  std::shared_ptr<PSP> psp;
  if (spl) psp = std::make_shared<PSPspl>(argv[optind], verbose);
  else     psp = std::make_shared<PSPout>(argv[optind], verbose);

  psp->PrintSummary(cout, stats, timeonly);

  return 0;
}
  
