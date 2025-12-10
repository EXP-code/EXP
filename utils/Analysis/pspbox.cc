/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in PSP files for a run and compute particles inside/outside
 *  spherical or cylindrical boxes
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

using namespace std;

                                // System libs
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include "numerical.H"
#include "ParticleReader.H"
#include "interp.H"
#include "massmodel.H"
#include "SphSL.H"
				// Support headers
#include "localmpi.H"
#include "cxxopts.H"
#include "foarray.H"
#include "libvars.H"


enum ProjectionType {Cylindrical=1, Spherical=2};

int main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, RMAX, ZMIN, ZMAX;
  int RBINS, ZBINS, IBEG, IEND, PBEG, PEND, ISKIP;
  std::string OUTFILE, INFILE, COMP, PROJ;
  std::string fileType, filePrefix;
  std::string psfiles, delim;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  cxxopts::Options options(argv[0], "Compute disk potential, force and density profiles\nfrom PSP phase-space output files\n");

  options.add_options()
    ("h,help", "Print this help message")
    ("v,verbose", "Print verbose processing info")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("psfiles", "Phase-space file list",
     cxxopts::value<std::string>(psfiles))
    ("RMIN", "minimum radius for output",
     cxxopts::value<double>(RMIN)->default_value("0.0"))
    ("RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("0.1"))
    ("ZMIN", "minimum z position",
     cxxopts::value<double>(ZMIN)->default_value("-1.0"))
    ("ZMAX", "maximum z position",
     cxxopts::value<double>(ZMAX)->default_value("1.0"))
    ("RBINS", "number of radial bins",
     cxxopts::value<int>(RBINS)->default_value("50"))
    ("ZBINS", "number of vertical bins",
     cxxopts::value<int>(ZBINS)->default_value("50"))
    ("PBEG", "first particle index",
     cxxopts::value<int>(PBEG)->default_value("0"))
    ("PEND", "last particle index",
     cxxopts::value<int>(PEND)->default_value("-1"))
    ("OUTFILE", "filename prefix",
     cxxopts::value<string>(OUTFILE)->default_value("gashisto"))
    ("COMP", "component name",
     cxxopts::value<string>(COMP)->default_value("dark"))
    ("PROJ", "Projection (cylindrical or spherical)",
     cxxopts::value<string>(PROJ)->default_value("cylindrical"))
    ;

  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }


  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    if (myid==0) std::cout << std::endl << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }

  bool verbose = false;
  if (vm.count("verbose")) verbose = true;


  // ==================================================
  // Get projection
  // ==================================================

  ProjectionType proj = Spherical;
  if (vm.count("PROJ")) {
    if (PROJ.compare("Cylindrical") == 0) proj = Cylindrical;
  }

  // ==================================================
  // Do round robin grid assignment of nodes
  // ==================================================

  auto files = PR::ParticleReader::parseFileList(psfiles, delim);
  
  unsigned nfiles = files.size();

  vector<double>   times0  (nfiles);
  vector<unsigned> inside0 (nfiles, 0), outside0 (nfiles, 0);
  vector<double>   minside0(nfiles, 0), moutside0(nfiles, 0);

  vector< vector<double> > linside0 (nfiles, vector<double>(3, 0));
  vector< vector<double> > loutside0(nfiles, vector<double>(3, 0));

  vector<double>   times  (nfiles);
  vector<unsigned> inside (nfiles, 0), outside (nfiles, 0);
  vector<double>   minside(nfiles, 0), moutside(nfiles, 0);

  vector< vector<double> > linside (nfiles, vector<double>(3, 0));
  vector< vector<double> > loutside(nfiles, vector<double>(3, 0));

  double ZCENTER = 0.5*(ZMAX + ZMIN);
  double ZWIDTH  = ZMAX - ZMIN;

  unsigned n = 0;

  for (auto batch : files) {

    PR::PRptr reader = PR::ParticleReader::createReader
      (fileType, batch, myid, true);

    times[n] = reader->CurrentTime();

    if (verbose and myid==0)
      std::cout << "Reading batch with T=" << times[n] << std::endl;

    // Find the component
    reader->SelectType(COMP);

    vector<double> L(3);
    int icnt = 0;
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {
      
      if (icnt > PBEG) {
	
	if (proj==Cylindrical) {
	  if (p->pos[2] >= ZCENTER-ZWIDTH && p->pos[2] <= ZCENTER+ZWIDTH) {
	    double R = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	    if (R>RMIN && R<RMAX) {
	      inside[n]++;
	      minside[n] += p->mass;
	      
	      L[0] = p->mass*
		(p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	      L[1] = p->mass*
		(p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
	      L[2] = p->mass*
		(p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);
	      
	      for (int k=0; k<3; k++) linside[n][k] += L[k];
	      
	    } else {
	      outside[n]++;
	      moutside[n] += p->mass;
	      
	      L[0] = p->mass*
		(p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	      L[1] = p->mass*
		(p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
	      L[2] = p->mass*
		(p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);
	      
	      for (int k=0; k<3; k++) loutside[n][k] += L[k];
	    }
	  }
	}
	else { // proj==Spherical
	  double R = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	  if (R>RMIN && R<RMAX) {
	    inside [n]++;
	    minside[n] += p->mass;
	    L[0] = p->mass*
	      (p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	    L[1] = p->mass*
	      (p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
	    L[2] = p->mass*
	      (p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);
	    
	    for (int k=0; k<3; k++) linside[n][k] += L[k];
	  } else {
	    outside [n]++;
	    moutside[n] += p->mass;
	    L[0] = p->mass*
	      (p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	    L[1] = p->mass*
	      (p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
	    L[2] = p->mass*
	      (p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);
	    
	    for (int k=0; k<3; k++) loutside[n][k] += L[k];
	  }
	}
      }
      
      if (PEND>0 && icnt>PEND) break;
      p = reader->nextParticle();
      icnt++;
    }
    if (verbose and myid==0)
      std::cout << "Processed " << icnt << " particles from batch " << n
		<< std::endl;

    n++;
  }
  
  MPI_Reduce(&times[0], &times0[0], nfiles,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(&inside[0], &inside0[0], nfiles,
	     MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&outside[0], &outside0[0], nfiles,
	     MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&minside[0], &minside0[0], nfiles,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&moutside[0], &moutside0[0], nfiles,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  for (int n=0; n<nfiles; n++) {
    MPI_Reduce(&linside[n][0], &linside0[n][0], 3,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&loutside[n][0], &loutside0[n][0], 3,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (myid==0) {
    string outf = OUTFILE + ".inout.counts";
    ofstream out(outf.c_str());
    
    out << left << setw(18) << "# Time"
	<< setw(12) << "| count in"
	<< setw(12) << "| count out"
	<< setw(18) << "| mass in"
	<< setw(18) << "| mass out"
	<< setw(18) << "| x angmom in"
	<< setw(18) << "| y angmom in"
	<< setw(18) << "| z angmom in"
	<< setw(18) << "| x angmom out"
	<< setw(18) << "| y angmom out"
	<< setw(18) << "| z angmom out"
	<< setw(18) << "| x angmom/M in"
	<< setw(18) << "| y angmom/M in"
	<< setw(18) << "| z angmom/M in"
	<< setw(18) << "| x angmom/M out"
	<< setw(18) << "| y angmom/M out"
	<< setw(18) << "| z angmom/M out"
	<< endl << setw(18) << "# 1"
	<< setw(12) << "| 2"
	<< setw(12) << "| 3"
	<< setw(18) << "| 4"
	<< setw(18) << "| 5"
	<< setw(18) << "| 6"
	<< setw(18) << "| 7"
	<< setw(18) << "| 8"
	<< setw(18) << "| 9"
	<< setw(18) << "| 10"
	<< setw(18) << "| 11"
	<< setw(18) << "| 12"
	<< setw(18) << "| 13"
	<< setw(18) << "| 14"
	<< setw(18) << "| 15"
	<< setw(18) << "| 16"
	<< setw(18) << "| 17"
	<< endl << right;
    
    for (int n=0; n<nfiles; n++) {
      out << setw(18) << times0   [n] 
	  << setw(12) << inside0  [n]
	  << setw(12) << outside0 [n]
	  << setw(18) << minside0 [n]
	  << setw(18) << moutside0[n];
      for (int k=0; k<3; k++)
	out << setw(18) << linside0[n][k];
      for (int k=0; k<3; k++)
	out << setw(18) << loutside0[n][k];
      for (int k=0; k<3; k++) {
	if (minside0[n]>0) out << setw(18) << linside0[n][k]/minside0[n];
	else out << setw(18) << 0.0;
      }
      for (int k=0; k<3; k++) {
	if (moutside0[n]>0) out << setw(18) << loutside0[n][k]/moutside0[n];
	else out << setw(18) << 0.0;
      }
      out << endl;
    }
    
  }

  MPI_Finalize();

  return 0;
}

