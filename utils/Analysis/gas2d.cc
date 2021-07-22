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
  
  int IBEG, IEND, NBINS, PBEG, PEND, ISKIP;
  double RMAX, ZCENTER, ZWIDTH;
  bool LOG;
  std::string OUTFILE, INFILE, RUNTAG, cname;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  po::options_description desc("Compute disk potential, force and density profiles\nfrom PSP phase-space output files\n\nAllowed options");
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("RMAX",                po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("ZCENTER",             po::value<double>(&ZCENTER)->default_value(0.0),
     "gas disk midplane")
    ("ZWIDTH",              po::value<double>(&ZWIDTH)->default_value(0.05),
     "gas disk halfwidth")
    ("NBINS",               po::value<int>(&NBINS)->default_value(100),
     "number of bins")
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
    ("LOG",                 po::value<bool>(&LOG)->default_value(false),
     "use logarithmic scaling for radial axis")
    ("OUTFILE",             po::value<string>(&OUTFILE)->default_value("gashisto"),
     "filename prefix")
    ("CNAME",               po::value<string>(&cname)->default_value("gas"),
     "Component name")
    ("INFILE",              po::value<string>(&INFILE)->default_value("OUT"),
     "phase space file")
    ("RUNTAG",              po::value<string>(&RUNTAG)->default_value("run"),
     "file containing desired indices for PSP output")
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

  double dR = 2.0*RMAX/(NBINS-1);

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
	  histo[n][k] = vector<double>(NBINS*NBINS, 0.0);

	// Do we need to close and reopen?
	if (in.rdstate() & ios::eofbit) {
	  in.close();
	  in.open(files[n].c_str());
	}

	int icnt = 0;
	vector<Particle> particles;

	PSPstanza *gas = psp->GetNamed(cname);
	if (!gas) {
	  if (myid==0) std::cerr << "No component named <" << cname << ">"
				 << std::endl;
	  MPI_Finalize();
	  exit(-2);
	}
	SParticle *p = psp->GetParticle();
	
	while (p) {

	  if (icnt > PBEG) {
	    if (p->pos(2) >= ZCENTER-ZWIDTH && p->pos(2) <= ZCENTER+ZWIDTH) {
	      int indX = static_cast<int>(floor( (p->pos(0) + RMAX)/dR ));
	      int indY = static_cast<int>(floor( (p->pos(1) + RMAX)/dR ));
	      if (indX >=0 && indX<NBINS &&
		  indY >=0 && indY<NBINS ) {
		histo[n][0][indY*NBINS+indX] += p->mass();
		histo[n][1][indY*NBINS+indX] += p->mass() * p->datr(0);
		histo[n][2][indY*NBINS+indX] += p->mass() * p->datr(1);
		histo[n][3][indY*NBINS+indX] += 
		  p->mass() * p->datr(0) * p->datr(1);
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
	histo[n][k] = vector<double>(NBINS*NBINS);
      for (int k=0; k<nval; k++)
	MPI_Recv(&histo[n][k][0], NBINS*NBINS, MPI_DOUBLE, 
		 curid, 12+k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (curid==myid) {
      double time;
      MPI_Send(&times[n], 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
      for (int k=0; k<nval; k++)
	MPI_Send(&histo[n][k][0], NBINS*NBINS, MPI_DOUBLE, 
		 0, 12+k, MPI_COMM_WORLD);
    }
  }
    
  if (myid==0) {
    ofstream out(OUTFILE + ".dat");

    out.precision(8);
    
    out << setw(8) << nval << setw(8) << NBINS << setw(8) << times.size()
	<< setw(18) << RMAX << endl;
    for (int n=0; n<nfiles; n++) {
      if (times.find(n)==times.end()) continue;
      out << left << setw(18) << times[n];
    }
    out << endl;
    for (int j=0; j<NBINS; j++)
      out << left << setw(18) << -RMAX + (0.5+j)*dR;
    out << endl;
    for (int n=0; n<nfiles; n++) {
      if (times.find(n)==times.end()) continue;
      for (int i=0; i<nval; i++) {
	for (unsigned j=0; j<NBINS; j++) {
	  for (unsigned k=0; k<NBINS; k++)
	    out << setw(18) << histo[n][i][j*NBINS+k];
	  out << endl;
	}
      }
    }
  }

  MPI_Finalize();

  return 0;
}

