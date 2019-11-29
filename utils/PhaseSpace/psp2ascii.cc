/*
  Separate a psp structure to ascii components

  MDWeinberg 06/10/02, 11/24/19
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
  cerr << prog << ": [-t time -v -h] filename\n\n";
  cerr << "    -o name         prefix name for each component (default: comp)\n";
  cerr << "    -a              use input format\n";
  cerr << "    -d dir          replacement SPL file directory\n";
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
  bool input   = false;
  string cname("comp"), new_dir("./");

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:o:avh");

    if (c == -1) break;

    switch (c) {

    case 'v':
      verbose = true;
      break;

    case 'o':
      cname.erase();
      cname = string(optarg);
      break;

    case 'd':
      new_dir.erase();
      new_dir = string(optarg);
      break;

    case 'a':
      input = true;
      break;

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

  std::string filename;

  if (optind < argc) {

    filename = std::string(argv[optind]);
    
    std::ifstream in(filename);
    if (!in) {
      cerr << "Error opening file <" << filename << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << filename << endl;

  } else {
    Usage(prog);
  }

				// Parse the PSP file
				// ------------------
  PSPptr psp;
  if (filename.find("SPL") != std::string::npos)
    psp = std::make_shared<PSPspl>(filename, new_dir);
  else
    psp = std::make_shared<PSPout>(filename);

				// Now write a summary
				// -------------------
  if (verbose) {

    psp->PrintSummary(cerr);
    
    cerr << "\nPSP file <" << filename << "> has time <" 
	 << psp->CurrentTime() << ">\n";
  }

				// Dump ascii for each component
				// -----------------------------
  PSPstanza *stanza;
  SParticle* part;

  for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
				// Open an output file for this stanza
				// -----------------------------------
    ostringstream oname;
    oname << cname << "." << stanza->name << '\0';
    ofstream out(oname.str().c_str());
    out.setf(ios::scientific);
    out.precision(10);

    if (!out) {
      cerr << "Couldn't open output name <" << oname.str() << ">\n";
      exit(-1);
    }
				// Print the header
				// ----------------
    out << setw(15) << stanza->comp.nbod 
	<< setw(10) << stanza->comp.niatr 
	<< setw(10) << stanza->comp.ndatr 
	<< endl;

    for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

      if (stanza->index_size) out << std::setw(18) << part->indx();
      out << std::setw(18) << part->mass();
      for (int i=0; i<3; i++) out << std::setw(18) << part->pos(i);
      for (int i=0; i<3; i++) out << std::setw(18) << part->vel(i);
      if (not input) out << std::setw(18) << part->phi();
      for (int i=0; i<part->niatr(); i++) out << std::setw(12) << part->iatr(i);
      for (int i=0; i<part->ndatr(); i++) out << std::setw(18) << part->datr(i);

      out << std::endl;		// End the record
    }
    
  }
  
  return 0;
}
  
