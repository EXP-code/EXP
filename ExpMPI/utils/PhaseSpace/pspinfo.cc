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

#include <PSP.H>

				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-v -t] filename\n";
  cerr << "        -v    verbose output\n";
  cerr << "        -t    print tipsy info\n";
  cerr << "        -T    print time info only\n";
  exit(-1);
}

int
main(int argc, char *argv[])
{
  bool tipsy = false;
  bool stats = false;
  bool timeonly = false;
  bool verbose = false;
  int c;
  
  while (1) {
    c = getopt(argc, argv, "vstTh");
    if (c == -1) break;

    switch (c) {
      
    case 'v':
      verbose = true;
      break;

    case 't':
      tipsy = true;
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

  cerr << "Filename: " << argv[optind] << endl;
  ifstream* in = new ifstream(argv[optind]);

  PSPDump psp(in, tipsy, verbose);
  psp.PrintSummary(in, cout, stats, timeonly);

  return 0;
}
  
