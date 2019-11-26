/*
  Separate a psp structure and make a histogram

  MDWeinberg 03/15/10, 11/24/19
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
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -x xmin         minimum x component\n";
  cerr << "    -X xmax         maximum x component\n";
  cerr << "    -x ymin         minimum y component\n";
  cerr << "    -Y ymax         maximum y component\n";
  cerr << "    -z zmin         minimum z component\n";
  cerr << "    -Z zmax         maximum z component\n";
  cerr << "    -1 numx         number of bins in x direction\n";
  cerr << "    -2 numy         number of bins in y direction\n";
  cerr << "    -o name         component name (default: comp)\n";
  cerr << "    -c comp         value index\n";
  cerr << "    -d dir          replacement SPL file directory\n";
  cerr << "    -r              spherical radius\n";
  cerr << "    -R              cylindrical radius\n";
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
  std::string cname("comp"), new_dir("");
  double xmin = 0.0, xmax = 1.0;
  double ymin = 0.0, ymax = 1.0;
  double zmin = -100.0, zmax = 100.0;
  int numx = 40;
  int numy = 40;
  int comp = 9;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "x:X:y:Y:z:Z:1:2:c:o:mnAvrRh");

    if (c == -1) break;

    switch (c) {

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

    case 'd':
      new_dir.erase();
      new_dir = string(optarg);
      break;

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

  std::string file;

  if (optind < argc) {

    file = std::string(argv[optind]);
    std::ifstream in(file);
    if (!in) {
      std::cerr << "Error opening file <" << file << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << file << endl;

  } else {
    Usage(prog);
  }


				// Parse the PSP file
				// ------------------
  PSPptr psp;
  if (file.find("SPL") != std::string::npos)
    psp = std::make_shared<PSPspl>(file, new_dir);
  else
    psp = std::make_shared<PSPout>(file);


				// Now write a summary
				// -------------------
  if (verbose) {

    psp->PrintSummary(cerr);
    
    cerr << "\nPSP file named <" << file << "> has time <" 
	 << psp->CurrentTime() << ">\n";
  }

				// Dump ascii for each component
				// -----------------------------
  
  double rtmp, mass, fac, val, dx=(xmax - xmin)/numx, dy=(ymax - ymin)/numy;
  vector<double> pos(3), vel(3);
  int itmp, ix, iy;

				// Make the array
				// --------------

  vector< vector<float> > value(numy), bmass(numy);
  for (int j=0; j<numy; j++) {
    value[j] = vector<float>(numx, 0);
    bmass[j] = vector<float>(numx, 0);
  }

  PSPstanza *stanza;
  SParticle* part;

  for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
    if (stanza->name != cname) continue;

    for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

      if (part->pos(0)<xmin || part->pos(0)>=xmax) continue;
      if (part->pos(1)<ymin || part->pos(1)>=ymax) continue;
      if (part->pos(2)<zmin || part->pos(2)>=zmax) continue;

      ix = static_cast<int>( floor( (part->pos(0) - xmin)/dx ) );
      iy = static_cast<int>( floor( (part->pos(1) - ymin)/dy ) );
      
      if (mweight) {
	bmass[iy][ix] += part->mass();
	fac = mass;
      } else {
	bmass[iy][ix] += 1.0;
	fac = 1.0;
      }

      if (comp == 0)
	value[iy][ix] += fac*part->mass();
      else if (comp <= 3)
	value[iy][ix] += fac*part->pos(comp-1);
      else if (comp <= 6)
	value[iy][ix] += fac*part->vel(comp-4);
      else if (comp <= 7 + part->niatr())
	value[iy][ix] += fac*part->iatr(comp-7);
      else
	value[iy][ix] += fac*part->datr(comp-7-part->niatr());
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
  
