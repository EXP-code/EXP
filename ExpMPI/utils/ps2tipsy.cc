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

#include <StringTok.H>
#include <header.H>

#define MAXDIM 3

typedef float Real;

struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
} ;

struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
} ;

struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
} ;

struct tipsydump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

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
  
  int ngas, ndark, nstar, ntot;
  list<Stanza> gas, dark, star;

  Dump () : ngas(0), ndark(0), nstar(0), ntot(0) {}
};

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool verbose = false;

  // Parse command line

  int c;
  int digit_optind = 0;

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


    for (int i=0; i<dump.header.ncomp; i++) {

      Stanza stanza;
      stanza.pos = in->tellg();
      
      ComponentHeader headerC;
      if (!headerC.read(in)) {
	cerr << "Error reading header\n";
	exit(-1);
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

				// Count up Tipsy types and make
				// linked  lists
				// -----------------------------
      if (!stanza.ttype.compare("gas")) {
	dump.ngas += stanza.nbod;
	dump.ntot += stanza.nbod;
	dump.gas.push_back(stanza);
      }
      if (!stanza.ttype.compare("dark")) {
	dump.ndark += stanza.nbod;
	dump.ntot += stanza.nbod;
	dump.dark.push_back(stanza);
      }
      if (!stanza.ttype.compare("star")) {
	dump.nstar += stanza.nbod;
	dump.ntot += stanza.nbod;
	dump.star.push_back(stanza);
      }

    }

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
      cerr << "          Gas particles:  " << itd->ngas << endl;
      cerr << "         Dark particles:  " << itd->ndark << endl;
      cerr << "         Star particles:  " << itd->nstar << endl;

      int cnt=1;

      for (its = itd->stanzas.begin(); its != itd->stanzas.end(); its++) {
	
				// Print the info for this stanza
				// ------------------------------
	cerr << setw(60) << setfill('-') << "-" << endl << setfill(' ');
	cerr << "--- Component #" << setw(2) << cnt++ << endl;
	cerr << setw(20) << " name :: "  << its->name   << endl
	     << setw(20) << " id :: "    << its->id     << endl
	     << setw(20) << " param :: " << its->param  << endl
	     << setw(20) << " tipsy :: " << its->ttype  << endl
	     << setw(20) << " nbod :: "  << its->nbod  << endl
	     << setw(20) << " niatr :: " << its->niatr << endl
	     << setw(20) << " ndatr :: " << its->ndatr << endl;
	cerr << setw(60) << setfill('-') << "-" << endl << setfill(' ');
	
      }
    }
  }
    
  cerr << "\nBest fit dump to <" << time << "> has time <" 
       << fid.header.time << ">\n";

				// Create tipsy output
				// -------------------
  in->close();
  delete in;
  in = new ifstream(argv[optind]);

  {
    tipsydump theader;
    gas_particle gas;
    dark_particle dark;
    star_particle star;

    theader.time = fid.header.time;
    theader.nbodies = fid.ntot;
    theader.ndim = 3;
    theader.nsph = fid.ngas;
    theader.ndark = fid.ndark;
    theader.nstar = fid.nstar;

    cout.write(&theader, sizeof(tipsydump));
    
    double rtmp;
    int itmp;

    list<Stanza>::iterator its;

				// Do gas particles
				// ----------------

    for (its = fid.gas.begin(); its != fid.gas.end(); its++) {

				// Position to beginning of particles
      in->seekg(its->pspos);

      for (int i=0; i<its->nbod; i++) {
	in->read(&rtmp, sizeof(double));
	gas.mass = rtmp;
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  gas.pos[i] = rtmp;
	}
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  gas.vel[i] = rtmp;
	}
	in->read(&rtmp, sizeof(double));
	gas.phi = rtmp;
	for (int i=0; i<its->niatr; i++) {
	  in->read(&itmp, sizeof(double));
	}
	for (int i=0; i<its->ndatr; i++) {
	  in->read(&rtmp, sizeof(double));
	}      

	gas.rho = gas.temp = gas.hsmooth = gas.metals = 0.0;

	cout.write(&gas, sizeof(gas_particle));
      }
    }

				// Do dark particles
				// -----------------

    for (its = fid.dark.begin(); its != fid.dark.end(); its++) {

				// Position to beginning of particles
      in->seekg(its->pspos);

      for (int i=0; i<its->nbod; i++) {
	in->read(&rtmp, sizeof(double));
	dark.mass = rtmp;
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  dark.pos[i] = rtmp;
	}
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  dark.vel[i] = rtmp;
	}
	in->read(&rtmp, sizeof(double));
	dark.phi = rtmp;
	for (int i=0; i<its->niatr; i++) {
	  in->read(&itmp, sizeof(double));
	}
	for (int i=0; i<its->ndatr; i++) {
	  in->read(&rtmp, sizeof(double));
	}      

	dark.eps = 0.0;

	cout.write(&dark, sizeof(dark_particle));

      }

    }

				// Do star particles
				// -----------------

    for (its = fid.star.begin(); its != fid.star.end(); its++) {

				// Position to beginning of particles
      in->seekg(its->pspos);

      for (int i=0; i<its->nbod; i++) {
	in->read(&rtmp, sizeof(double));
	star.mass = rtmp;
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  star.pos[i] = rtmp;
	}
	for (int i=0; i<3; i++) {
	  in->read(&rtmp, sizeof(double));
	  star.vel[i] = rtmp;
	}
	in->read(&rtmp, sizeof(double));
	star.phi = rtmp;
	for (int i=0; i<its->niatr; i++) {
	  in->read(&itmp, sizeof(double));
	}
	for (int i=0; i<its->ndatr; i++) {
	  in->read(&rtmp, sizeof(double));
	}      

	star.metals = star.tform = star.eps = 0.0;

	cout.write(&star, sizeof(star_particle));
      }
    }
    
  }
  
  return 0;
}
  
