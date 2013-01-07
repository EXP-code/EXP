/*
  Separate a psp structure and make a histogram

  MDWeinberg 03/15/10
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
  cerr << "    -x xmin         minimum x component\n";
  cerr << "    -X xmax         maximum x component\n";
  cerr << "    -x ymin         minimum y component\n";
  cerr << "    -Y ymax         maximum y component\n";
  cerr << "    -z zmin         minimum z component\n";
  cerr << "    -Z zmax         maximum z component\n";
  cerr << "    -1 numx         number of bins in x direction\n";
  cerr << "    -2 numy         number of bins in y direction\n";
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
  double xmin = 0.0, xmax = 1.0;
  double ymin = 0.0, ymax = 1.0;
  double zmin = -100.0, zmax = 100.0;
  int numx = 40;
  int numy = 40;
  int comp = 9;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:x:X:y:Y:z:Z:1:2:c:o:mnAvh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 'x':
      xmin = atof(optarg);
      break;

    case 'X':
      xmax = atof(optarg);
      break;

    case 'y':
      ymin = atof(optarg);
      break;

    case 'Y':
      ymax = atof(optarg);
      break;

    case 'z':
      zmin = atof(optarg);
      break;

    case 'Z':
      zmax = atof(optarg);
      break;

    case '1':
      numx = atoi(optarg);
      break;

    case '2':
      numy = atoi(optarg);
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
  double rtmp, mass, fac, val, dx=(xmax - xmin)/numx, dy=(ymax - ymin)/numy;
  vector<double> pos(3), vel(3);
  int itmp, icnt, ix, iy;

				// Make the array
				// --------------

  vector< vector<float> > value(numy), bmass(numy);
  for (int j=0; j<numy; j++) {
    value[j] = vector<float>(numx, 0);
    bmass[j] = vector<float>(numx, 0);
  }

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

      if (pos[0]<xmin || pos[0]>=xmax) continue;
      if (pos[1]<ymin || pos[1]>=ymax) continue;
      if (pos[2]<zmin || pos[2]>=zmax) continue;

      ix = static_cast<int>( floor( (pos[0] - xmin)/dx ) );
      iy = static_cast<int>( floor( (pos[1] - ymin)/dy ) );
      
      if (mweight) {
	bmass[iy][ix] += mass;
	fac = mass;
      } else {
	bmass[iy][ix] += 1.0;
	fac = 1.0;
      }

      if (comp == 0)
	value[iy][ix] += fac*mass;
      else if (comp <= 3)
	value[iy][ix] += fac*pos[comp-1];
      else if (comp <= 6)
	value[iy][ix] += fac*vel[comp-4];
      else
	value[iy][ix] += fac*vals[comp-7];
    }
    
  }
  
  //
  // Output
  //
  cout.write((const char *)&numx, sizeof(int));
  cout.write((const char *)&numy, sizeof(int));

  float f, m=0.0, rhomin=1e20, rhomax=0.0;

  for (int i=0; i<numx; i++) {
    f = xmin + dx*(0.5+i);
    cout.write((const char *)&f, sizeof(float));
  }

  for (int j=0; j<numy; j++) {
    f = ymin + dy*(0.5+j);
    cout.write((const char *)&f, sizeof(float));
  }

  for (int j=0; j<numy; j++) {
    for (int i=0; i<numx; i++) {
      m += bmass[j][i];
      if (areal)  {
	f = value[j][i]/dx*dy;
      } else {
	if (bmass[j][i] > 0.0) f = value[j][i]/bmass[j][i];
	else f = 0.0;
      }
      cout.write((const char *)&f, sizeof(float));
      rhomin = min<float>(f, rhomin);
      rhomax = max<float>(f, rhomax);
    }
  }

  cerr << "Total mass=" << m 
       << ", minD=" << rhomin 
       << ", maxD=" << rhomax << endl;

  return 0;
}
  
