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
string outdir, runtag;

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-t time -4 -8 -s -v -h] filename\n\n";
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -s              add a cparam string (prompts user)\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -4              use float for real values\n";
  cerr << "    -8              use double for real values (default)\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}

int
main(int argc, char **argv)
{
  char *prog      = argv[0];
  double time     = 1e20;
  bool use_float  = false;
  bool add_cparam = false;
  bool verbose    = false;
  

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:s48vh");

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

    case '4':
      use_float = true;
      break;

    case '8':
      use_float = false;
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
    
    PSPstanza *stanza;
    SParticle* part;

    for (stanza=psp.GetStanza(); stanza!=0; stanza=psp.NextStanza()) {

				// Position to header
      in->seekg(stanza->pos);
      
      if (add_cparam) {

	cerr <<"===================================================" << endl
	     << "Name=" << stanza->name << endl
	     << "ID=" << stanza->id << endl
	     << "Current cparam string=" << stanza->cparam << endl
	     << "Enter new cparam string: ";
	char line[1024];
	cin.getline(line, 1024);
	stanza->cparam = line;

	string delim = " : ";
	string infostr = 
	  stanza->name + delim + 
	  stanza->id + delim + 
	  stanza->cparam + delim + 
	  stanza->fparam + '\0';
	if (infostr.size() > stanza->comp.ninfochar) {
	  stanza->comp.ninfochar = infostr.size() + 1;
	  stanza->comp.info = 
	    boost::shared_array<char>(new char [stanza->comp.ninfochar]);
	}
	strcpy(stanza->comp.info.get(), infostr.c_str());
      }

      stanza->comp.write(&std::cout);

				// Position to beginning of particles
      in->seekg(stanza->pspos);

				// Write the particles
      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in))
	part->write(&std::cout, use_float, stanza->index_size);
    }

  }
  
  return 0;
}
  
