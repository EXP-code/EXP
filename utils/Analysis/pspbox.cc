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
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>
				// Support headers
#include <localmpi.H>
#include <cxxopts.H>
#include <foarray.H>
#include <libvars.H>


enum ProjectionType {Cylindrical=1, Spherical=2};

int main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, RMAX, ZMIN, ZMAX;
  int RBINS, ZBINS, IBEG, IEND, PBEG, PEND, ISKIP;
  std::string OUTFILE, INFILE, RUNTAG, COMP, PROJ;
  std::string fileType, filePrefix, runtag, dir;
  bool GNUPLOT;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  cxxopts::Options options(argv[0], "Compute disk potential, force and density profiles\nfrom PSP phase-space output files\n");

  options.add_options()
    ("h,help", "Print this help message")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("P,prefix", "prefix for phase-space files",
     cxxopts::value<std::string>(filePrefix)->default_value("OUT"))
    ("D,dir", "directory for phase-space files",
     cxxopts::value<std::string>(dir)->default_value("."))
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
    ("IBEG", "first PSP index",
     cxxopts::value<int>(IBEG)->default_value("0"))
    ("IEND", "last PSP index",
     cxxopts::value<int>(IEND)->default_value("100"))
    ("ISKIP", "skip PSP interval",
     cxxopts::value<int>(ISKIP)->default_value("1"))
    ("PBEG", "first particle index",
     cxxopts::value<int>(PBEG)->default_value("0"))
    ("PEND", "last particle index",
     cxxopts::value<int>(PEND)->default_value("-1"))
    ("OUTFILE", "filename prefix",
     cxxopts::value<string>(OUTFILE)->default_value("gashisto"))
    ("INFILE", "phase space file",
     cxxopts::value<string>(INFILE)->default_value("OUT"))
    ("RUNTAG", "file containing desired indices for PSP output",
     cxxopts::value<string>(RUNTAG)->default_value("run"))
    ("COMP", "component name",
     cxxopts::value<string>(COMP)->default_value("dark"))
    ("PROJ", "Projection (cylindrical or spherical)",
     cxxopts::value<string>(PROJ)->default_value("cylindrical"))
    ("GNUPLOT", "Write gnuplot type output",
     cxxopts::value<bool>(GNUPLOT)->default_value("false"))
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

  std::vector<string> files;
				// Root nodes looks for existence of files
  if (myid==0) {
    for (int i=IBEG; i<=IEND; i++) {
      auto lab = PR::ParticleReader::fileNameCreator
	(fileType, i, myid, dir, runtag, filePrefix);
      std::ifstream in(lab);
      if (in) files.push_back(lab);
      else break;
      std::cout << "." << i << flush;
    }
    std::cout << endl;
  }

  unsigned nfiles = files.size();
  MPI_Bcast(&nfiles, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  if (nfiles==0) {
    if (myid==0) std::cout << "No files found!" << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  for (unsigned n=0; n<nfiles; n++) {
    unsigned sz;
    if (myid==0) sz = files[n].size();
    MPI_Bcast(&sz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (myid==0) {
      char *c = const_cast<char*>(files[n].c_str());
      MPI_Bcast(c, sz+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
      char l[sz+1];
      MPI_Bcast(&l[0], sz+1, MPI_CHAR, 0, MPI_COMM_WORLD);
      files.push_back(l);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

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

  for (int n=0; n<nfiles; n++) {

    PR::PRptr reader = PR::ParticleReader::createReader
      (fileType, files[n], myid, true);

    times[n] = reader->CurrentTime();

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
    string outf = RUNTAG + ".inout.counts";
    ofstream out(outf.c_str());
    
    out << left << setw(18) << "# Time"
	<< setw(18) << "| count in"
	<< setw(18) << "| count out"
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
	<< setw(18) << "| 2"
	<< setw(18) << "| 3"
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
	  << setw(10) << inside0  [n]
	  << setw(10) << outside0 [n]
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

