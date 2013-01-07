/*
  Separate a psp structure and make a 1-d histogram

  MDWeinberg 08/26/11
*/

using namespace std;

#include <cstdlib>

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
string outdir, runtag;

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-t time -v -h] filename\n\n";
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -a axis         x=1, y=2, z=3 (default: 3)\n";
  cerr << "    -p pmin         minimum position along axis\n";
  cerr << "    -P pmax         maximum position along axis\n";
  cerr << "    -b numb         number of bins\n";
  cerr << "    -c comp         value index\n";
  cerr << "    -o name         component name (default: comp)\n";
  cerr << "    -A              areal average\n";
  cerr << "    -m              mass-weighted values\n";
  cerr << "    -n              number-weighted values\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool mweight = true;
  bool nweight = false;
  bool areal   = false;
  bool verbose = false;
  string cname("comp");
  double pmin = -100.0, pmax = 100.0;
  int axis = 3;
  int numb = 40;
  int comp = 9;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:a:p:P:b:c:o:mnAvh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 'a':
      axis = atoi(optarg);
      break;

    case 'p':
      pmin = atof(optarg);
      break;

    case 'P':
      pmax = atof(optarg);
      break;

    case 'b':
      numb = atoi(optarg);
      break;

    case 'c':
      comp = atoi(optarg);
      break;

    case 'm':
      mweight = true;
      nweight = false;
      break;

    case 'n':
      mweight = false;
      nweight = true;
      break;

    case 'A':
      areal = true;
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

  } else {
    Usage(prog);
  }

				// Axis sanity check 
				// ------------------
  if (axis<1) axis = 1;
  if (axis>3) axis = 3;

				// Parse the PSP file
				// ------------------
  PSPDump psp(in);

  in->close();

				// Now write a summary
				// -------------------
  if (verbose) {

    psp.PrintSummary(in, cerr);
    
    cerr << "\nBest fit dump to <" << time << "> has time <" 
	 << psp.SetTime(time) << ">\n";
  } else 
    psp.SetTime(time);

				// Dump ascii for each component
				// -----------------------------
  delete in;
  in = new ifstream(argv[optind]);

  
  list<PSPstanza>::iterator its;
  double rtmp, mass, fac, dp=(pmax - pmin)/numb;
  vector<double> pos(3), vel(3);
  int itmp, icnt, iv;

				// Make the array
				// --------------

  vector<float> value(numb, 0), bmass(numb, 0);

  for (its = psp.CurrentDump()->stanzas.begin(); 
       its != psp.CurrentDump()->stanzas.end(); its++) {

    if (its->name != cname) continue;


				// Position to beginning of particles
    in->seekg(its->pspos);

    icnt = 0;
    vector<double> vals;

    for (int j=0; j<its->comp.nbod; j++) {
      in->read((char *)&mass, sizeof(double));
      for (int i=0; i<3; i++) {
	  in->read((char *)&pos[i], sizeof(double));
      }
      for (int i=0; i<3; i++) {
	  in->read((char *)&vel[i], sizeof(double));
      }
      in->read((char *)&rtmp, sizeof(double));
      vals.push_back(rtmp);
      for (int i=0; i<its->comp.niatr; i++) {
	in->read((char *)&itmp, sizeof(int));
	vals.push_back(itmp);
      }
      for (int i=0; i<its->comp.ndatr; i++) {
	in->read((char *)&rtmp, sizeof(double));
	vals.push_back(rtmp);
      }      

      if (pos[axis-1]<pmin || pos[axis-1]>=pmax) continue;

      iv = static_cast<int>( floor( (pos[axis-1] - pmin)/dp ) );
      
      if (mweight) {
	bmass[iv] += mass;
	fac = mass;
      } else {
	bmass[iv] += 1.0;
	fac = 1.0;
      }

      if (comp == 0)
	value[iv] += fac*mass;
      else if (comp <= 3)
	value[iv] += fac*pos[comp-1];
      else if (comp <= 6)
	value[iv] += fac*vel[comp-4];
      else
	value[iv] += fac*vals[comp-7];
    }
    
  }
  
  //
  // Output
  //
  double Time = psp.CurrentTime();

  float p, f, m=0.0;

  for (int i=0; i<numb; i++) {
    p  = pmin + dp*(0.5+i);
    f  = 0.0;
    m += bmass[i];
    if (areal)  {
      f = value[i]/dp;
    } else {
      if (bmass[i] > 0.0) f = value[i]/bmass[i];
    }
    cout << setw(18) << Time 
	 << setw(18) << p
	 << setw(18) << f
	 << setw(18) << m
	 << endl;
  }
  cout << endl;

  return 0;
}
