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

				// Boost stuff
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace po = boost::program_options;


                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <PSP.H>

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
std::string outdir, outfile, runtag;
boost::mt19937 random_gen;

//=============================================================================

enum ComponentType {Star=1, Gas=2, Halo=4};

int    CFLAGS;
int    TEMP;
int    DENS;
int    KNUD;
int    STRL;
int    Nbins;
int    AXIS;
double Rmin;
double Rmax;

void add_particles(PSPptr psp, vector< vector<double> >& ret, int& nhist)
{
  
  int nbods = psp->GetStanza()->comp.nbod;

  nhist = 4;

  if (CFLAGS & Gas) {
    if (TEMP>=0) nhist++;
    if (DENS>=0) nhist++;
    if (KNUD>=0) nhist++;
    if (STRL>=0) nhist++;
  }

  ret = vector< vector<double> >(Nbins);

  SParticle *part = psp->GetParticle();
  double pos;
  int indx;
    
  for (int i=0; i<nbods; i++) {
    if (part==0) {
      cerr << "Error reading particle [n=" << 0 << ", i=" << i << "]" << endl;
      exit(-1);
    }

    pos = part->pos(AXIS);

    if (pos>=Rmin && pos<Rmax) {
    
      indx = static_cast<int>(floor((pos-Rmin)*Nbins/(Rmax - Rmin)));
      if (ret[indx].size() == 0) ret[indx] = vector<double>(nhist, 0.0);

      int cnt = 0;

      ret[indx][cnt++] += part->mass();
      ret[indx][cnt++] += part->mass() * part->vel(0);
      ret[indx][cnt++] += part->mass() * part->vel(1);
      ret[indx][cnt++] += part->mass() * part->vel(2);

      if (CFLAGS & Gas) {
	if (TEMP>=0) ret[indx][cnt++] += part->mass() * part->datr(TEMP);
	if (DENS>=0) ret[indx][cnt++] += part->mass() * part->datr(DENS);
	if (KNUD>=0) ret[indx][cnt++] += part->mass() * part->datr(KNUD);
	if (STRL>=0) ret[indx][cnt++] += part->mass() * part->datr(STRL);
      }
    }

    part = psp->NextParticle();
  }

  return;
}


int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  int IMIN, IMAX;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  po::options_description desc("Compute 1-dimensional projection of shocktube runs\nAllowed options");
  desc.add_options()
    ("help,h",                                                                       "Print this help message")
    ("OUT",
     "assume that PSP files are in original format")
    ("SPL",
     "assume that PSP files are in split format")
    ("CFLAGS",              po::value<int>(&CFLAGS)->default_value(2),
     "component flags (Star=1)")
    ("TEMP",                po::value<int>(&TEMP)->default_value(0),
     "temperature (default=0)")
    ("DENS",                po::value<int>(&DENS)->default_value(1),
     "density (default=1)")
    ("KNUD",                po::value<int>(&KNUD)->default_value(4),
     "Knudsen (default=4)")
    ("STRL",                po::value<int>(&STRL)->default_value(5),
     "Straoul (default=5)")
    ("Rmin",                po::value<double>(&Rmin)->default_value(0.0),
     "minimum position")
    ("Rmax",                po::value<double>(&Rmax)->default_value(1.0),
     "maximum position")
    ("Nbins",               po::value<int>(&Nbins)->default_value(100),
     "number of bins")
    ("AXIS",                po::value<int>(&AXIS)->default_value(2),
     "which axis")
    ("OUTFILE",             po::value<string>(&outfile)->default_value("slab.prof"),
     "output filename")
    ("RUNTAG",              po::value<string>(&runtag)->default_value("run"),
     "run tag")
    ("MININDX",             po::value<int>(&IMIN)->default_value(0),
     "Minimum PSP index")
    ("MAXINDX",             po::value<int>(&IMAX)->default_value(100),
     "Maximum PSP index")
    ;
  

  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }


  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    std::cout << std::endl << desc << std::endl;
    return 0;
  }

  // ==================================================
  // Open frame list
  // ==================================================

  ofstream indx;
  ifstream in;

  ofstream out(outfile);
  if (!out) {
    cerr << "Error opening output file <" << outfile
	 << ">" << endl;
    exit(-1);
  }

  for (int i=IMIN; i<=IMAX; i++) {

    ostringstream sout;
    if (vm.count("SPL")) sout << "SPL.";
    else                 sout << "OUT.";
    sout << runtag << "." << right << setw(5) << setfill('0') << i;

    in.close();
    in.open(sout.str().c_str());
    if (!in) continue;

    cout << "Reading <" << sout.str() << "> . . ." << flush;
  
    PSPptr psp;
    if (vm.count("SPL")) psp = std::make_shared<PSPspl>(sout.str());
    else                 psp = std::make_shared<PSPout>(sout.str());

    double time = psp->CurrentTime();

    vector< vector<double> > ret;
    int nhist;
    double dz = (Rmax-Rmin)/Nbins;

    add_particles(psp, ret, nhist);
      
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
	
    cout << " done" << endl;
  }

  return 0;
}

