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

				// EXP classes
#include <ParticleReader.H>
#include <cxxopts.H>
#include <libvars.H>

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

void add_particles(PR::PRptr reader, std::vector< std::vector<double> >& ret,
		   int& nhist)
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
  
  std::string fileType, filePrefix, outfile, runtag;
  std::string psfiles, delim;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  cxxopts::Options options(argv[0], "Compute 1-dimensional projection of shocktube runs\n");

  options.add_options()
    ("h,help", "Print this help message")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("P,prefix", "prefix for phase-space files",
     cxxopts::value<std::string>(filePrefix)->default_value("OUT"))
    ("compname", "train on Component (default=slab)",
     cxxopts::value<std::string>(cname)->default_value("slab"))
    ("TEMP", "temperature (default=0)",
     cxxopts::value<int>(TEMP)->default_value("0"))
    ("DENS", "density (default=1)",
     cxxopts::value<int>(DENS)->default_value("1"))
    ("KNUD", "Knudsen (default=4)",
     cxxopts::value<int>(KNUD)->default_value("4"))
    ("STRL", "Straoul (default=5)",
     cxxopts::value<int>(STRL)->default_value("5"))
    ("Rmin", "minimum position",
     cxxopts::value<double>(Rmin)->default_value("0.0"))
    ("Rmax", "maximum position",
     cxxopts::value<double>(Rmax)->default_value("1.0"))
    ("Nbins", "number of bins",
     cxxopts::value<int>(Nbins)->default_value("100"))
    ("AXIS", "which axis",
     cxxopts::value<int>(AXIS)->default_value("2"))
    ("OUTFILE", "output filename",
     cxxopts::value<string>(outfile)->default_value("slab.prof"))
    ("RUNTAG", "run tag",
     cxxopts::value<string>(runtag)->default_value("run"))
    ("psfile", "List of phase space files for processing",
     cxxopts::value<std::string>(psfiles))
    ("delimiter", "Phase-space file list delimiter for node index",
     cxxopts::value<std::string>(delim))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    std::cout << std::endl << options.help() << std::endl;
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

  auto files = PR::ParticleReader::parseFileList(psfiles, delim);

  for (auto batch : files) {

    PR::PRptr reader = PR::ParticleReader::createReader
      (fileType, batch, myid, true);

    reader->SelectType(cname);

    double time = reader->CurrentTime();

    std::vector< std::vector<double> > ret;
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

