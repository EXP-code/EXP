// This is really -*- C++ -*-

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Open PSP dumps and compute one-dimensional traces (histograms)
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/28/08
 *
 ***************************************************************************/

				// C++/STL headers
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <PSP.H>
#include <ProgramParam.H>

//=============================================================================
// Variables not used but needed for linking
//=============================================================================
//
int VERBOSE = 4;
int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
double tpos = 0.0;
double tnow = 0.0;
int myid = 0;  
//
//=============================================================================
// Option database
//=============================================================================

program_option init[] = {
  {"CFLAGS",		"int",		"2",		"component flags (Star=1, Gas=2, Halo=4)"},
  {"TEMP",		"int",		"0",		"temperature (default=0)"},
  {"DENS",		"int",		"1",		"density (default=1)"},
  {"Rmin",		"double",	"0.0",		"minimum position"},
  {"Rmax",		"double",	"1.0",		"maximum position"},
  {"Nbins",		"int",		"100",		"number of bins"},
  {"AXIS",		"int",		"2",		"which axis"},
  {"OUTFILE",		"string",	"slab.prof",	"output filename"},
  {"RUNTAG",		"string",	"run",		"run tag"},
  {"MININDX",		"int",		"0",		"Minimum PSP index"},
  {"MAXINDX",		"int",		"100",		"Maximum PSP index"},
  {"",			"",		"",		""}
};


const char desc[] = "Compute 1-dimensional projection of shocktube runs\n";

ProgramParam config(desc, init);

//=============================================================================

enum ComponentType {Star=1, Gas=2, Halo=4};

int    CFLAGS;
int    TEMP;
int    DENS;
int    Nbins;
int    AXIS;
double Rmin;
double Rmax;

void add_particles(ifstream* in, PSPDump* psp,
		   vector< vector<double> >& ret, int& nhist)
{
  
  int nbods = 0;

  nhist = 4;

  if (CFLAGS & Star) {
    nbods = psp->CurrentDump()->nstar;
    psp->GetStar();
  }

  if (CFLAGS & Gas) {
    nbods = psp->CurrentDump()->ngas;
    psp->GetGas();
    if (TEMP>=0) nhist++;
    if (DENS>=0) nhist++;
  }

  if (CFLAGS & Halo) {
    nbods = psp->CurrentDump()->ndark;
    psp->GetDark();
  }

  ret = vector< vector<double> >(Nbins);

  SParticle *part = psp->GetParticle(in);
  double pos;
  int indx;
    
  for (int i=0; i<nbods; i++) {
    if (part==0) {
      cerr << "Error reading particle [n=" << 0 << ", i=" << i << "]" << endl;
      exit(-1);
    }

    pos = part->pos[AXIS];

    if (pos>=Rmin && pos<Rmax) {
    
      indx = static_cast<int>(floor((pos-Rmin)*Nbins/(Rmax - Rmin)));
      if (ret[indx].size() == 0) ret[indx] = vector<double>(nhist, 0.0);

      int cnt = 0;

      ret[indx][cnt++] += part->mass;
      ret[indx][cnt++] += part->mass * part->vel[0];
      ret[indx][cnt++] += part->mass * part->vel[1];
      ret[indx][cnt++] += part->mass * part->vel[2];

      if (CFLAGS & Gas) {
	if (TEMP>=0) ret[indx][cnt++] += part->mass * part->datr[TEMP];
	if (DENS>=0) ret[indx][cnt++] += part->mass * part->datr[DENS];
      }
    }

    part = psp->NextParticle(in);
  }

  return;
}


int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  if (config.parse_args(argc,argv)) return -1;
  
  
  CFLAGS  = config.get<int>("CFLAGS");
  TEMP    = config.get<int>("TEMP");
  DENS    = config.get<int>("DENS");
  Nbins   = config.get<int>("Nbins");
  AXIS    = config.get<int>("AXIS");
  Rmin    = config.get<double>("Rmin");
  Rmax    = config.get<double>("Rmax");

  // ==================================================
  // Open frame list
  // ==================================================

  ofstream indx;
  ifstream in;

  int IMIN = config.get<int>("MININDX");
  int IMAX = config.get<int>("MAXINDX");

  ofstream out(config.get<string>("OUTFILE").c_str());
  if (!out) {
    cerr << "Error opening output file <" << config.get<string>("OUTFILE")
	 << ">" << endl;
    exit(-1);
  }

  for (int i=IMIN; i<=IMAX; i++) {

    ostringstream sout;
    sout << "OUT." << config.get<string>("RUNTAG") << "."
	 << right << setw(5) << setfill('0') << i;

    in.close();
    in.open(sout.str().c_str());
    if (!in) continue;

    cout << "Reading <" << sout.str() << "> . . ." << flush;
  
    PSPDump psp(&in, true);

    Dump *dump = psp.GetDump();
    
    while (dump) {

      double time = dump->header.time;

      if (in.rdstate() & ios::eofbit) {
	in.close();
	in.open(sout.str().c_str());
      }

      vector< vector<double> > ret;
      int nhist;
      double dz = (Rmax-Rmin)/Nbins;

      add_particles(&in, &psp, ret, nhist);
      
      for (int n=0; n<Nbins; n++) {
	out << setw(15) << time
	    << setw(15) << Rmin + dz*(0.5+n);
	if (ret[n].size()) {
	  out << setw(15) << ret[n][0]/dz;
	  for (int j=1; j<nhist; j++) out << setw(15) << ret[n][j]/ret[n][0];
	}
	else
	  for (int j=0; j<nhist; j++) out << setw(15) << 0.0;

	out << endl;
      }
      out << endl;
	
      dump = psp.NextDump();
    }

    cout << " done" << endl;
  }

  return 0;
}

