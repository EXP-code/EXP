/*
  Info on each stanza of a phase space dump

  MDWeinberg 03/17/02
*/

#include <unistd.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>

#include <PSP.H>

#define MAXDIM 3


				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-t time -v -h] filename\n\n";
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool verbose = false;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:vh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
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

  ifstream *in;

  if (optind < argc) {

    ifstream *in2 = new ifstream(argv[optind]);
    if (!*in2) {
      cerr << "Error opening file <" << argv[optind] << "> for input\n";
      exit(-1);
    }

    cerr << "Using filename: " << argv[optind] << endl;

				// Assign file stream to input stream
    in = in2;

  }

  PSPDump psp(in, true);
  in->close();
  delete in;

				// Now write a summary
				// -------------------

  if (verbose) psp.PrintSummary(cerr);

    
  cerr << "\nBest fit dump to <" << time << "> has time <" 
       << psp.SetTime(time) << ">\n";

				// Create new dump
				// ---------------
  in = new ifstream(argv[optind]);

  {
    cout.write(&psp.CurrentDump()->header, sizeof(MasterHeader));
    
    double rtmp;
    int itmp;

    list<Stanza>::iterator its;

    for (its = psp.CurrentDump()->stanzas.begin(); 
	 its != psp.CurrentDump()->stanzas.end(); its++) {

				// Position to header
      in->seekg(its->pos);
      
      ComponentHeader headerC;
      if (!headerC.read(in)) {
	cerr << "Error reading header\n";
	exit(-1);
      }
      headerC.write(&cout);

				// Position to beginning of particles
      in->seekg(its->pspos);

      for (int i=0; i<its->nbod; i++) {
	in->read(&rtmp, sizeof(double));
	cout.write(&rtmp, sizeof(double));
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  cout.write(&rtmp, sizeof(double));
	}
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  cout.write(&rtmp, sizeof(double));
	}
	in->read(&rtmp, sizeof(double));
	cout.write(&rtmp, sizeof(double));
	for (int i=0; i<its->niatr; i++) {
	  in->read(&itmp, sizeof(double));
	  cout.write(&itmp, sizeof(int));
	}
	for (int i=0; i<its->ndatr; i++) {
	  in->read(&rtmp, sizeof(double));
	  cout.write(&rtmp, sizeof(double));
	}      
      }
    }

  }
  
  return 0;
}
  
