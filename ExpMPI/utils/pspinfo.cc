/*
  Info on each stanza of a phase space dump

  MDWeinberg 03/17/02
*/

#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

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
  cerr << prog << ": filename\n";
  exit(-1);
}

int
main(int argv, char *argc[])
{
  streampos lastpos, currpos;
  vector<streampos> poslist;

  MasterHeader header;

  if (argv == 1) Usage(argc[0]);

  cerr << "Filename: " << argc[1] << endl;
  ifstream* in = new ifstream(argc[1]);

  while (1) {
    currpos = in->tellg();

				// Read the header, quit on failure
				// --------------------------------
    if(!in->read(&header, sizeof(MasterHeader))) break;
    lastpos = currpos;
    poslist.push_back(currpos);

    cout << "Time=" << header.time << "   [" << currpos << "]" << endl;
    cout << "   Total particle number: " << header.ntot  << endl;
    cout << "   Number of components:  " << header.ncomp << endl;

    for (int i=0; i<header.ncomp; i++) {

      ComponentHeader headerC;

      if (!headerC.read(in)) {
	cerr << "Error reading header\n";
	exit(-1);
      }

				// Parse the info string
				// ---------------------
      StringTok<string> tokens(headerC.info);
      string name = trimLeft(trimRight(tokens(":")));
      string id = trimLeft(trimRight(tokens(":")));
      string param = trimLeft(trimRight(tokens(":")));

				// Strip of the tipsy type
      StringTok<string> tipsytype(name);
      string ttype = trimLeft(trimRight(tipsytype(" ")));

				// Print the info for this stanza
				// ------------------------------
      cout << setw(60) << setfill('-') << "-" << endl << setfill(' ');
      cout << "--- Component #" << setw(2) << i+1 << endl;
      cout << setw(20) << " name :: "  << name          << endl
	   << setw(20) << " id :: "    << id            << endl
	   << setw(20) << " param :: " << param         << endl
	   << setw(20) << " tipsy :: " << ttype         << endl
	   << setw(20) << " nbod :: "  << headerC.nbod  << endl
	   << setw(20) << " niatr :: " << headerC.niatr << endl
	   << setw(20) << " ndatr :: " << headerC.ndatr << endl;
      cout << setw(60) << setfill('-') << "-" << endl << setfill(' ');
      
				// Skip forward to next header
				// ---------------------------
      in->seekg(headerC.nbod*(8*sizeof(double)             + 
			      headerC.niatr*sizeof(int)    +
			      headerC.ndatr*sizeof(double)
			      ), ios::cur);
    }

  }

}
  
