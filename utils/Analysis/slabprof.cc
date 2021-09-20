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

namespace po = boost::program_options;


                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// EXP classes
#include <ParticleReader.H>

//=============================================================================
// Global variables
//=============================================================================

int         TEMP;
int         DENS;
int         KNUD;
int         STRL;
int         Nbins;
int         AXIS;
double      Rmin;
double      Rmax;
std::string cname;

void add_particles(PRptr reader, vector< vector<double> >& ret, int& nhist)
{
  
  nhist = 4;

  if (cname.find("Gas") == 0) {
    if (TEMP>=0) nhist++;
    if (DENS>=0) nhist++;
    if (KNUD>=0) nhist++;
    if (STRL>=0) nhist++;
  }

  ret = vector< vector<double> >(Nbins);

  double pos;
  int indx;
    
  for (auto part=reader->firstParticle(); part!=0; part=reader->nextParticle()) {

    pos = part->pos[AXIS];

    if (pos>=Rmin && pos<Rmax) {
    
      indx = static_cast<int>(floor((pos-Rmin)*Nbins/(Rmax - Rmin)));
      if (ret[indx].size() == 0) ret[indx] = vector<double>(nhist, 0.0);

      int cnt = 0;

      ret[indx][cnt++] += part->mass;
      ret[indx][cnt++] += part->mass * part->vel[0];
      ret[indx][cnt++] += part->mass * part->vel[1];
      ret[indx][cnt++] += part->mass * part->vel[2];

      if (cname.find("Gas") == 0) {
	if (TEMP>=0) ret[indx][cnt++] += part->mass * part->dattrib[TEMP];
	if (DENS>=0) ret[indx][cnt++] += part->mass * part->dattrib[DENS];
	if (KNUD>=0) ret[indx][cnt++] += part->mass * part->dattrib[KNUD];
	if (STRL>=0) ret[indx][cnt++] += part->mass * part->dattrib[STRL];
      }
    }
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
  std::string fileType, filePrefix, outfile, runtag;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  po::options_description desc("Compute 1-dimensional projection of shocktube runs\nAllowed options");
  desc.add_options()
    ("help,h",
     "Print this help message")
    ("filetype,F",
     po::value<std::string>(&fileType)->default_value("PSPout"),
     "input file type")
    ("prefix,P",
     po::value<std::string>(&filePrefix)->default_value("OUT"),
     "prefix for phase-space files")
    ("compname",            po::value<std::string>(&cname)->default_value("slab"),
     "train on Component (default=slab)")
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

    auto file = ParticleReader::fileNameCreator(fileType, i, "", runtag);
    
    PRptr reader = ParticleReader::createReader(fileType, file, true);
    reader->SelectType(cname);

    double time = reader->CurrentTime();

    vector< vector<double> > ret;
    int nhist;
    double dz = (Rmax-Rmin)/Nbins;

    add_particles(reader, ret, nhist);
      
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

