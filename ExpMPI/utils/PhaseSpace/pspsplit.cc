/*
  Info on each stanza of a phase space dump

  MDWeinberg 03/17/02
*/

#include <cstdlib>
#include <cstring>
#include <cmath>
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
  cerr << prog << ": [-t time -s -v -h] filename\n\n";
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -s              add a cparam string (prompts user)\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool add_cparam = false;
  bool verbose = false;
  

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:svh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 's':
      add_cparam = true;
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

  if (verbose) psp.PrintSummary(in, cerr);

    
  cerr << "\nBest fit dump to <" << time << "> has time <" 
       << psp.SetTime(time) << ">\n";

				// Create new dump
				// ---------------
  in = new ifstream(argv[optind]);

  {
    cout.write((char *)&psp.CurrentDump()->header, sizeof(MasterHeader));
    
    double rtmp;
    int itmp;

    list<PSPstanza>::iterator its;

    for (its = psp.CurrentDump()->stanzas.begin(); 
	 its != psp.CurrentDump()->stanzas.end(); its++) {

				// Position to header
      in->seekg(its->pos);
      
      ComponentHeader headerC;
      if (!headerC.read(in)) {
	cerr << "Error reading header\n";
	exit(-1);
      }

      if (add_cparam) {

	cerr <<"===================================================" << endl
	     << "Name=" << its->name << endl
	     << "ID=" << its->id << endl
	     << "Current cparam string=" << its->cparam << endl
	     << "Enter new cparam string: ";
	char line[1024];
	cin.getline(line, 1024);
	its->cparam = line;

	string delim = " : ";
	string infostr = 
	  its->name + delim + 
	  its->id + delim + 
	  its->cparam + delim + 
	  its->fparam + '\0';
	if (infostr.size() > headerC.ninfochar) {
	  delete [] headerC.info;
	  headerC.ninfochar = infostr.size() + 1;
	  headerC.info = new char [headerC.ninfochar];
	}
	strcpy(headerC.info, infostr.c_str());
      }

      headerC.write(&cout);

				// Position to beginning of particles
      in->seekg(its->pspos);

      for (int i=0; i<its->nbod; i++) {
	in->read((char *)&rtmp, sizeof(double));
	cout.write((char *)&rtmp, sizeof(double));
	for (int i=0; i<3; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  cout.write((char *)&rtmp, sizeof(double));
	}
	for (int i=0; i<3; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  cout.write((char *)&rtmp, sizeof(double));
	}
	in->read((char *)&rtmp, sizeof(double));
	cout.write((char *)&rtmp, sizeof(double));
	for (int i=0; i<its->niatr; i++) {
	  in->read((char *)&itmp, sizeof(double));
	  cout.write((char *)&itmp, sizeof(int));
	}
	for (int i=0; i<its->ndatr; i++) {
	  in->read((char *)&rtmp, sizeof(double));
	  cout.write((char *)&rtmp, sizeof(double));
	}      
      }
    }

  }
  
  return 0;
}
  
