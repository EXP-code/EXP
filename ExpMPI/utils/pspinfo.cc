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
  cerr << prog << ": filename [tipsy info (true/false)]\n";
  exit(-1);
}

int
main(int argv, char *argc[])
{
  if (argv < 2 || argv > 3) Usage(argc[0]);

  cerr << "Filename: " << argc[1] << endl;
  ifstream* in = new ifstream(argc[1]);

  bool tipsy = false;
  if (argv == 3 && atoi(argc[2])) tipsy = true;

  PSPDump psp(in, tipsy);
  psp.PrintSummary(cout);

}
  
