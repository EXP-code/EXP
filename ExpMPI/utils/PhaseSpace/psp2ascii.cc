/*
  Separate a psp structure to ascii components

  MDWeinberg 06/10/02
*/

using namespace std;

#include <unistd.h>
#include <stdlib.h>

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

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-t time -v -h] filename\n\n";
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -o name         prefix name for each component (default: comp)\n";
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
  string cname("comp");

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:o:vh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 'v':
      verbose = true;
      break;

    case 'o':
      cname.erase();
      cname = string(optarg);
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

    if (verbose) cerr << "Using filename: " << argv[optind] << endl;

				// Assign file stream to input stream
    in = in2;

  }


				// Parse the PSP file
				// ------------------
  PSPDump psp(in);

  in->close();
  delete in;

				// Now write a summary
				// -------------------
  if (verbose) {

    psp.PrintSummary(cerr);
    
    cerr << "\nBest fit dump to <" << time << "> has time <" 
	 << psp.SetTime(time) << ">\n";
  } else 
    psp.SetTime(time);

				// Dump ascii for each component
				// -----------------------------
  in = new ifstream(argv[optind]);

  
  list<PSPstanza>::iterator its;
  double rtmp;
  int itmp;

  for (its = psp.CurrentDump()->stanzas.begin(); 
       its != psp.CurrentDump()->stanzas.end(); its++) {

				// Open an output file
				// -------------------

    ostringstream oname;
    oname << cname << "." << its->name << '\0';
    ofstream out(oname.str().c_str());
    out.setf(ios::scientific);
    out.precision(10);

    if (!out) {
      cerr << "Couldn't open output name <" << oname.str() << ">\n";
      exit(-1);
    }

				// Print the header

    out << setw(15) << its->nbod 
	<< setw(10) << its->niatr 
	<< setw(10) << its->ndatr 
	<< endl;

				// Position to beginning of particles
    in->seekg(its->pspos);

    for (int j=0; j<its->nbod; j++) {
      in->read((char *)&rtmp, sizeof(double));
      out << setw(18) << rtmp;
      for (int i=0; i<3; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  out << setw(18) << rtmp;
      }
      for (int i=0; i<3; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  out << setw(18) << rtmp;
      }
      in->read((char *)&rtmp, sizeof(double));
      out << setw(18) << rtmp;
      for (int i=0; i<its->niatr; i++) {
	in->read((char *)&itmp, sizeof(int));
	out << setw(12) << itmp;
      }
      for (int i=0; i<its->ndatr; i++) {
	in->read((char *)&rtmp, sizeof(double));
	out << setw(12) << rtmp;
      }      

      out << endl;		// End the record

    }
    
  }
  
  return 0;
}
  
