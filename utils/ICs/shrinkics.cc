#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include <getopt.h>
#include <pthread.h>

// EXP library support
//
#include "libvars.H"

using namespace std;

//===========================================================================

void usage(char *prog)
{
  cout << "Routine to select the first N particles from an ascii PSP file" 
       << endl << endl << "Usage:" << endl << endl
       << prog << " [options]" << endl << endl
       << setw(15) << "Option" << setw(10) << "Argument" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Description" << endl
       << resetiosflags(ios::left)
       << endl
       << setw(15) << "-n or --number" << setw(10) << "0" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "New number of particles (if 0, numer is /10)" << endl
       << resetiosflags(ios::left)
       << resetiosflags(ios::left)
       << setw(15) << "-i or --file" << setw(10) << "halo.bods" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Body file" << endl
       << resetiosflags(ios::left)
       << endl;

  exit(0);
}

int main(int argc, char** argv)
{
  const int bufsiz = 2048;
  char buffer[bufsiz];
  
  int newnum = 0;
  string filename = "halo.bods";

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"number", 0, 0, 0},
      {"file", 1, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "n:i:h",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 0:
      {
	string optname(long_options[option_index].name);

	if (!optname.compare("number")) {
	  newnum = atoi(optarg);
	} else if (!optname.compare("file")) {
	  filename = optarg;
	} else {
	  cout << "Option " << long_options[option_index].name;
	  if (optarg) cout << " with arg " << optarg;
	  cout << " is not defined " << endl;
	  exit(0);
	}
      }
      break;

    case 'n':
      newnum = atoi(optarg);
      break;

    case 'i':
      filename = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
    }

  }

  //===================
  // Check input file
  //===================
  ifstream in(filename.c_str());
  if (!in) {
    cerr << "Error opening input file <" << filename << ">" << endl;
    exit(-1);
  }

  //===================
  // Read info line
  //===================
  
  int orignum, ndatr, niatr;
  in.getline(buffer, bufsiz);
  {
    istringstream sin(buffer);
    sin >> orignum;
    sin >> ndatr;
    sin >> niatr;
  }

  if (newnum == 0) newnum = orignum/10;

  double mfac = (double)orignum/newnum;
  vector<double> ps(7);
  vector<double> dd(ndatr);
  vector<int>    ii(niatr);

  cout << setw(12) << newnum
       << setw(12) << ndatr
       << setw(12) << niatr
       << endl;

  cout.setf(ios::scientific);
  cout.precision(12);

  double mtot = 0.0;

  for (int i=0; i<newnum; i++) {
    in.getline(buffer, bufsiz);
    istringstream sin(buffer);

    for (int k=0; k<7; k++)     sin >> ps[k];
    for (int k=0; k<ndatr; k++) sin >> dd[k];
    for (int k=0; k<niatr; k++) sin >> ii[k];
    
    ps[0] *= mfac;

    for (int k=0; k<7; k++)     cout << setw(20) << ps[k];
    for (int k=0; k<ndatr; k++) cout << setw(20) << dd[k];
    for (int k=0; k<niatr; k++) cout << setw(20) << ii[k];
    cout << endl;

    mtot += ps[0];
  }

  cerr << "Total mass = " << mtot << endl;

  return 0;
}




