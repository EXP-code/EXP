/*
  Separate a psp structure to ascii components

  MDWeinberg 06/10/02
*/

#include <unistd.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>

#include <StringTok.H>
#include <header.H>

extern string trimLeft(const string);
extern string trimRight(const string);

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

struct Stanza {
  streampos pos, pspos;
  string name;
  string id;
  string param;
  string ttype;
  int nbod;
  int niatr;
  int ndatr;
};

class Dump 
{
public:

  streampos pos;
  MasterHeader header;
  list<Stanza> stanzas;
  
  Dump () {}
};

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

    cerr << "Using filename: " << argv[optind] << endl;

				// Assign file stream to input stream
    in = in2;

  }
				// Look for best fit time
  double tdif = 1.0e30;
  Dump fid;

  list<Dump> dumps;

  while (1) {

    Dump dump;

    dump.pos = in->tellg();
				// Read the header, quit on failure
				// --------------------------------
    if(!in->read((char *)&dump.header, sizeof(MasterHeader))) break;


    bool ok = true;

    for (int i=0; i<dump.header.ncomp; i++) {

      Stanza stanza;
      stanza.pos = in->tellg();
      
      ComponentHeader headerC;
      if (!headerC.read(in)) {
	cerr << "Error reading component header\n";
	ok = false;
	break;
      }

      stanza.pspos = in->tellg();

				// Parse the info string
				// ---------------------
      StringTok<string> tokens(headerC.info);
      stanza.name = trimLeft(trimRight(tokens(":")));
      stanza.id = trimLeft(trimRight(tokens(":")));
      stanza.param = trimLeft(trimRight(tokens(":")));

				// Strip of the tipsy type
      StringTok<string> tipsytype(stanza.name);
      stanza.ttype = trimLeft(trimRight(tipsytype(" ")));
      stanza.nbod = headerC.nbod;
      stanza.niatr = headerC.niatr;
      stanza.ndatr = headerC.ndatr;

      
				// Skip forward to next header
				// ---------------------------
      in->seekg(headerC.nbod*(8*sizeof(double)             + 
			      headerC.niatr*sizeof(int)    +
			      headerC.ndatr*sizeof(double)
			      ), ios::cur);

      dump.stanzas.push_back(stanza);

    }

    if (!ok) break;

    dumps.push_back(dump);
    if (fabs(time - dump.header.time) < tdif) {
      fid = dump;
      tdif = fabs(time-dump.header.time);
    }
  }


				// Now write a summary
				// -------------------
  if (verbose) {

    list<Dump>::iterator itd;
    list<Stanza>::iterator its;

    for (itd = dumps.begin(); itd != dumps.end(); itd++) {

      cerr << "Time=" << itd->header.time << "   [" << itd->pos << "]" << endl;
      cerr << "   Total particle number: " << itd->header.ntot  << endl;
      cerr << "   Number of components:  " << itd->header.ncomp << endl;

      int cnt=1;

      for (its = itd->stanzas.begin(); its != itd->stanzas.end(); its++) {
	
				// Print the info for this stanza
				// ------------------------------
	cerr << setw(60) << setfill('-') << "-" << endl << setfill(' ');
	cerr << "--- Component #" << setw(2) << cnt++ << endl;
	cerr << setw(20) << " name :: "  << its->name   << endl
	     << setw(20) << " id :: "    << its->id     << endl
	     << setw(20) << " param :: " << its->param  << endl
	     << setw(20) << " nbod :: "  << its->nbod  << endl
	     << setw(20) << " niatr :: " << its->niatr << endl
	     << setw(20) << " ndatr :: " << its->ndatr << endl;
	cerr << setw(60) << setfill('-') << "-" << endl << setfill(' ');
	
      }
    }
  }
    
  cerr << "\nBest fit dump to <" << time << "> has time <" 
       << fid.header.time << ">\n";

				// Dump ascii for each component
				// -----------------------------
  in->close();
  delete in;

  


  in = new ifstream(argv[optind]);

  
  list<Stanza>::iterator its;
  double rtmp;
  int itmp;

  for (its = fid.stanzas.begin(); its != fid.stanzas.end(); its++) {

				// Open an output file
				// -------------------

    ostrstream oname;
    oname << cname << "." << its->name << '\0';
    ofstream out(oname.str());
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

    for (int i=0; i<its->nbod; i++) {
      in->read(&rtmp, sizeof(double));
      out << setw(18) << rtmp;
      for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  out << setw(18) << rtmp;
      }
      for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  out << setw(18) << rtmp;
      }
      in->read(&rtmp, sizeof(double));
      out << setw(18) << rtmp;
      for (int i=0; i<its->niatr; i++) {
	in->read(&itmp, sizeof(double));
	out << setw(12) << itmp;
      }
      for (int i=0; i<its->ndatr; i++) {
	in->read(&rtmp, sizeof(double));
	out << setw(12) << rtmp;
      }      

      out << endl;		// End the record

    }
    
  }
  
  return 0;
}
  
