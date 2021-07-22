/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in PSP files for a run and compute 2-d gas distribution
 *  histogram
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
//
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
  
  double RMAX, ZMIN, ZMAX;
  int RBINS, IBEG, IEND, ISKIP, PBEG, PEND, ZBINS;
  std::string OUTFILE, INFILE, RUNTAG, CNAME;
  bool GNUPLOT;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  po::options_description desc("\nCompute disk potential, force and density profiles\nfrom PSP phase-space output files\n\nAllowed options");
  desc.add_options()
    ("help,h",                                                                       "Print this help message")
    ("OUT",
     "assume that PSP files are in original format")
    ("SPL",
     "assume that PSP files are in split format")
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
    ("CNAME",               po::value<string>(&CNAME)->default_value("gas"),
     "Name for gas component")
    ("RUNTAG",              po::value<string>(&RUNTAG)->default_value("run"),
     "file containing desired indices for PSP output")
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

  double dR = RMAX/RBINS;
  double dZ = (ZMAX - ZMIN)/ZBINS;

  const int nval = 4;
  map<int, vector< vector<double> > > histo;
  map<int, double> times;

  for (int n=0; n<nfiles; n++) {

    if (n % numprocs == myid) {

      PSPptr psp;
      if (vm.count("SPL")) psp = std::make_shared<PSPspl>(files[n]);
      else                 psp = std::make_shared<PSPout>(files[n]);

      if (psp) {

	times[n] = psp->CurrentTime();
	histo[n] = vector< vector<double> >(nval);
	for (int k=0; k<nval; k++) 
	  histo[n][k] = vector<double>(RBINS*ZBINS, 0.0);

	// Do we need to close and reopen?
	if (in.rdstate() & ios::eofbit) {
	  in.close();
	  in.open(files[n].c_str());
	}

	int icnt = 0;
	vector<Particle> particles;

	PSPstanza *gas = psp->GetNamed(CNAME);
	if (!gas) {
	  if (myid==0) std::cerr << "No component named <" << CNAME << ">"
				 << std::endl;
	  MPI_Finalize();
	  exit(-2);
	}
	SParticle *p = psp->GetParticle();
	
	while (p) {

	  if (icnt > PBEG) {
	    double R = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1) );
	    if (p->pos(2) >= ZMIN && p->pos(2) < ZMAX && R < RMAX) {
	      int indR = static_cast<int>(floor( R/dR ));
	      int indZ = static_cast<int>(floor( (p->pos(2) - ZMIN)/dZ ));
	      if (indR >=0 && indR<RBINS &&
		  indZ >=0 && indZ<ZBINS ) {
		histo[n][0][indZ*RBINS+indR] += p->mass();
		histo[n][1][indZ*RBINS+indR] += p->mass() * p->datr(0);
		histo[n][2][indZ*RBINS+indR] += p->mass() * p->datr(1);
		histo[n][3][indZ*RBINS+indR] += p->mass() * p->datr(0) * p->datr(1);
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
  
  cout << "[#" << myid << ", " << times.size() <<  ", " 
       << histo.size() << "]" << endl;

  for (int n=0; n<nfiles; n++) {
    int curid = n % numprocs;
    if (curid==0) continue;
    if (myid==0) {
      double time;
      MPI_Recv(&time, 1, MPI_DOUBLE, curid, 11, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      times[n] = time;
      histo[n] = vector< vector<double> >(nval);
      for (int k=0; k<nval; k++) 
	histo[n][k] = vector<double>(RBINS*ZBINS);
      for (int k=0; k<nval; k++)
	MPI_Recv(&histo[n][k][0], RBINS*ZBINS, MPI_DOUBLE, 
		 curid, 12+k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (curid==myid) {
      double time;
      MPI_Send(&times[n], 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
      for (int k=0; k<nval; k++)
	MPI_Send(&histo[n][k][0], RBINS*ZBINS, MPI_DOUBLE, 
		 0, 12+k, MPI_COMM_WORLD);
    }
  }
    
  if (myid==0) {

    for (int n=0; n<nfiles; n++) {
      if (times.find(n)==times.end()) continue;

      ostringstream sout ;
      sout << OUTFILE << "." << n;
      ofstream out(sout.str().c_str());
      
      out.precision(8);
    
      if (GNUPLOT) {
	out << "# Time=" << times[n] << endl;

	for (unsigned j=0; j<ZBINS; j++) {
	  for (unsigned k=0; k<RBINS; k++) {
	    out << setw(18) << dR*(0.5+k) << setw(18) << ZMIN + dZ*(0.5+j);
	    out << setw(18) << histo[n][0][j*RBINS+k];
	    if (histo[n][0][j*RBINS+k]>0.0) {
	      for (int i=1; i<nval; i++)
		out << setw(18) 
		    << histo[n][i][j*RBINS+k]/histo[n][0][j*RBINS+k];
	    } else {
	      for (int i=1; i<nval; i++) out << setw(18) << 0.0;
	    }
	    out << endl;
	  }
	  out << endl;
	}
      } else {
	out << setw(18) << times[n] << endl
	    << setw(10) << RBINS << setw(10) << ZBINS << endl;
	for (unsigned k=0; k<RBINS; k++) out << setw(18) << dR*(0.5+k);
	out << endl;
	for (unsigned j=0; j<ZBINS; j++) out << setw(18) << ZMIN + dZ*(0.5+j);
	out << endl;
	for (unsigned k=0; k<RBINS; k++) {
	  for (unsigned j=0; j<ZBINS; j++) {
	    for (unsigned k=0; k<RBINS; k++) {
	      out << setw(18) << histo[n][0][j*RBINS+k];
	      if (histo[n][0][j*RBINS+k]>0.0) {
		for (int i=1; i<nval; i++)
		  out << setw(18) 
		      << histo[n][i][j*RBINS+k]/histo[n][0][j*RBINS+k];
	      } else {
		for (int i=1; i<nval; i++) out << setw(18) << 0.0;
	      }
	      out << endl;
	    }
	  }
	}
      }
    }
  }

  MPI_Finalize();

  return 0;
}

