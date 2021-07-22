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

				// Boost stuff

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace po = boost::program_options;


                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include "Particle.h"
#include <PSP.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>

#include <localmpi.H>
#include <foarray.H>


enum ProjectionType {Cylindrical=1, Spherical=2};

				// Variables not used but needed for linking
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
string outdir, runtag;
double tpos = 0.0;
double tnow = 0.0;
boost::mt19937 random_gen;
  
int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMIN, RMAX, ZMIN, ZMAX;
  int RBINS, ZBINS, IBEG, IEND, PBEG, PEND, ISKIP;
  std::string OUTFILE, INFILE, RUNTAG, COMP, PROJ;
  bool GNUPLOT;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  po::options_description desc("Compute disk potential, force and density profiles\nfrom PSP phase-space output files\n\nAllowed options");
  desc.add_options()
    ("help,h",                                                                       "Print this help message")
    ("OUT",
     "assume that PSP files are in original format")
    ("SPL",
     "assume that PSP files are in split format")
    ("RMIN",                po::value<double>(&RMAX)->default_value(0.0),
     "minimum radius for output")
    ("RMAX",                po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("ZMIN",                po::value<double>(&ZMIN)->default_value(-1.0),
     "minimum z position")
    ("ZMAX",                po::value<double>(&ZMAX)->default_value(1.0),
     "maximum z position")
    ("RBINS",               po::value<int>(&RBINS)->default_value(50),
     "number of radial bins")
    ("ZBINS",               po::value<int>(&ZBINS)->default_value(50),
     "number of vertical bins")
    ("IBEG",                po::value<int>(&IBEG)->default_value(0),
     "first PSP index")
    ("IEND",                po::value<int>(&IEND)->default_value(100),
     "last PSP index")
    ("ISKIP",               po::value<int>(&ISKIP)->default_value(1),
     "skip PSP interval")
    ("PBEG",                po::value<int>(&PBEG)->default_value(0),
     "first particle index")
    ("PEND",                po::value<int>(&PEND)->default_value(-1),
     "last particle index")
    ("OUTFILE",             po::value<string>(&OUTFILE)->default_value("gashisto"),
     "filename prefix")
    ("INFILE",              po::value<string>(&INFILE)->default_value("OUT"),
     "phase space file")
    ("RUNTAG",              po::value<string>(&RUNTAG)->default_value("run"),
     "file containing desired indices for PSP output")
    ("COMP",                po::value<string>(&COMP)->default_value("dark"),
     "component name")
    ("PROJ",                po::value<string>(&PROJ)->default_value("cylindrical"),
     "Projection (cylindrical or spherical)")
    ("GNUPLOT",             po::value<bool>(&GNUPLOT)->default_value(false),
     "Write gnuplot type output")
    ;

  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
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
  // Get projection
  // ==================================================

  ProjectionType proj = Spherical;
  if (vm.count("PROJ")) {
    if (PROJ.compare("Cylindrical") == 0) proj = Cylindrical;
  }

  // ==================================================
  // Do round robin grid assignment of nodes
  // ==================================================

  ofstream indx;
  ifstream in;

  vector<string> files;
				// Root nodes looks for existence of files
  if (myid==0) {
    for (int i=IBEG; i<=IEND; i++) {
      ostringstream lab;
      lab << INFILE << "." 
	  << RUNTAG << "." 
	  << setw(5) << right << setfill('0') << i;
      ifstream in(lab.str().c_str());
      if (in) files.push_back(lab.str());
      else break;
      cout << "." << i << flush;
    }
    cout << endl;
  }

  unsigned nfiles = files.size();
  MPI_Bcast(&nfiles, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
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

    if (n % numprocs == myid) {

      PSPptr psp;
      if (vm.count("SPL")) psp = std::make_shared<PSPspl>(files[n]);
      else                 psp = std::make_shared<PSPout>(files[n]);

      if (psp) {

	times[n] = psp->CurrentTime();

	// Do we need to close and reopen?
	if (in.rdstate() & ios::eofbit) {
	  in.close();
	  in.open(files[n].c_str());
	}

	// Find the component
	PSPstanza *stanza;
	for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
	  if (stanza->name == COMP) break;
	}

	if (stanza==0) {
	  std::cout << "Could not find Component <" << COMP << "> at time = "
		    << psp->CurrentTime() << std::endl;
	} else {

	  in.seekg(stanza->pspos);

	  vector<double> L(3);
	  int icnt = 0;
	  for (SParticle* 
		 p=psp->GetParticle(); p!=0; p=psp->NextParticle()) {

	    if (icnt > PBEG) {

	      if (proj==Cylindrical) {
		if (p->pos(2) >= ZCENTER-ZWIDTH && p->pos(2) <= ZCENTER+ZWIDTH) {
		  double R = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
		  if (R>RMIN && R<RMAX) {
		    inside[n]++;
		    minside[n] += p->mass();

		    L[0] = p->mass()*
		      (p->pos(1)*p->vel(2) - p->pos(2)*p->vel(1));
		    L[1] = p->mass()*
		      (p->pos(2)*p->vel(0) - p->pos(0)*p->vel(2));
		    L[2] = p->mass()*
		      (p->pos(0)*p->vel(1) - p->pos(1)*p->vel(0));
		    
		    for (int k=0; k<3; k++) linside[n][k] += L[k];

		  } else {
		    outside[n]++;
		    moutside[n] += p->mass();

		    L[0] = p->mass()*
		      (p->pos(1)*p->vel(2) - p->pos(2)*p->vel(1));
		    L[1] = p->mass()*
		      (p->pos(2)*p->vel(0) - p->pos(0)*p->vel(2));
		    L[2] = p->mass()*
		      (p->pos(0)*p->vel(1) - p->pos(1)*p->vel(0));
		    
		    for (int k=0; k<3; k++) loutside[n][k] += L[k];
		  }
		}
	      }
	      else { // proj==Spherical
		double R = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
		if (R>RMIN && R<RMAX) {
		  inside [n]++;
		  minside[n] += p->mass();
		  L[0] = p->mass()*
		    (p->pos(1)*p->vel(2) - p->pos(2)*p->vel(1));
		  L[1] = p->mass()*
		    (p->pos(2)*p->vel(0) - p->pos(0)*p->vel(2));
		  L[2] = p->mass()*
		    (p->pos(0)*p->vel(1) - p->pos(1)*p->vel(0));
		  
		  for (int k=0; k<3; k++) linside[n][k] += L[k];
		} else {
		  outside [n]++;
		  moutside[n] += p->mass();
		  L[0] = p->mass()*
		    (p->pos(1)*p->vel(2) - p->pos(2)*p->vel(1));
		  L[1] = p->mass()*
		    (p->pos(2)*p->vel(0) - p->pos(0)*p->vel(2));
		  L[2] = p->mass()*
		    (p->pos(0)*p->vel(1) - p->pos(1)*p->vel(0));
		  
		  for (int k=0; k<3; k++) loutside[n][k] += L[k];
		}
	      }
	    }
	    
	    if (PEND>0 && icnt>PEND) break;
	    p = psp->NextParticle();
	    icnt++;
	  }
	}
      }
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

