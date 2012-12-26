/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in PSP files for a run and compute
 *  2-d gas distribution histogram
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
#include <Vector.h>
#include <numerical.h>
#include <Particle.h>
#include <PSP.H>
#include <interp.h>
#include <massmodel.h>
#include <SphereSL.H>

#include <localmpi.h>
#include <ProgramParam.H>
#include <foarray.H>

program_option init[] = {
  {"RMAX",		"double",	"0.1",		"maximum radius for output"},
  {"ZCENTER",		"double",	"0.0",		"gas disk midplane"},
  {"ZWIDTH",		"double",	"0.05",		"gas disk halfwidth"},
  {"NBINS",		"int",		"100",		"number of bins"},
  {"IBEG",		"int",		"0",		"first PSP index"},
  {"IEND",		"int",		"100",		"last PSP index"},
  {"ISKIP",		"int",		"1",		"skip PSP interval"},
  {"PBEG",		"int",		"0",		"first particle index"},
  {"PEND",		"int",		"-1",		"last particle index"},
  {"LOG",		"bool",		"false",	"use logarithmic scaling for radial axis"},
  {"OUTFILE",		"string",	"gashisto",	"filename prefix"},
  {"INFILE",		"string",	"OUT",		"phase space file"},
  {"RUNTAG",		"string",	"run",		"file containing desired indices for PSP output"},
  {"",			"",		"",		""}
};


const char desc[] = "Compute disk potential, force and density profiles from PSP phase-space output files\n";


ProgramParam config(desc, init);

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
  
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  // ==================================================
  // Do round robin grid assignment of nodes
  // ==================================================

  ofstream indx;
  ifstream in;

  vector<string> files;
				// Root nodes looks for existence of files
  if (myid==0) {
    for (int i=config.get<int>("IBEG"); i<=config.get<int>("IEND"); i++) {
      ostringstream lab;
      lab << config.get<string>("INFILE") << "." 
	  << config.get<string>("RUNTAG") << "." 
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

  double rmax = config.get<double>("RMAX");
  int   nbins = config.get<int>   ("NBINS");
  double zcen = config.get<double>("ZCENTER");
  double zwid = config.get<double>("ZWIDTH");
  int    pbeg = config.get<int>   ("PBEG");
  int    pend = config.get<int>   ("PEND");

  double dR = 2.0*rmax/(nbins-1);

  const int nval = 4;
  map<int, vector< vector<double> > > histo;
  map<int, double> times;

  for (int n=0; n<nfiles; n++) {

    if (n % numprocs == myid) {

      ifstream in(files[n].c_str());
      PSPDump psp(&in, true);

      Dump *dump = psp.GetDump();
      
      if (dump) {

	times[n] = psp.CurrentTime();
	histo[n] = vector< vector<double> >(nval);
	for (int k=0; k<nval; k++) 
	  histo[n][k] = vector<double>(nbins*nbins, 0.0);

	// Do we need to close and reopen?
	if (in.rdstate() & ios::eofbit) {
	  in.close();
	  in.open(files[n].c_str());
	}

	int icnt = 0;
	vector<Particle> particles;

	PSPstanza *gas = psp.GetGas();
	SParticle *p = psp.GetParticle(&in);
	
	while (p) {

	  if (icnt > pbeg) {
	    if (p->pos(2) >= zcen-zwid && p->pos(2) <= zcen+zwid) {
	      int indX = static_cast<int>(floor( (p->pos(0) + rmax)/dR ));
	      int indY = static_cast<int>(floor( (p->pos(1) + rmax)/dR ));
	      if (indX >=0 && indX<nbins &&
		  indY >=0 && indY<nbins ) {
		histo[n][0][indY*nbins+indX] += p->mass();
		histo[n][1][indY*nbins+indX] += p->mass() * p->datr(0);
		histo[n][2][indY*nbins+indX] += p->mass() * p->datr(1);
		histo[n][3][indY*nbins+indX] += 
		  p->mass() * p->datr(0) * p->datr(1);
	      }
	    }
	  }
	    
	  if (pend>0 && icnt>pend) break;
	  p = psp.NextParticle(&in);
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
	histo[n][k] = vector<double>(nbins*nbins);
      for (int k=0; k<nval; k++)
	MPI_Recv(&histo[n][k][0], nbins*nbins, MPI_DOUBLE, 
		 curid, 12+k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (curid==myid) {
      double time;
      MPI_Send(&times[n], 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
      for (int k=0; k<nval; k++)
	MPI_Send(&histo[n][k][0], nbins*nbins, MPI_DOUBLE, 
		 0, 12+k, MPI_COMM_WORLD);
    }
  }
    
  if (myid==0) {
    string outf = config.get<string>("OUTFILE") + ".dat";
    ofstream out(outf.c_str());

    out.precision(8);
    
    out << setw(8) << nval << setw(8) << nbins << setw(8) << times.size()
	<< setw(18) << rmax << endl;
    for (int n=0; n<nfiles; n++) {
      if (times.find(n)==times.end()) continue;
      out << left << setw(18) << times[n];
    }
    out << endl;
    for (int j=0; j<nbins; j++)
      out << left << setw(18) << -rmax + (0.5+j)*dR;
    out << endl;
    for (int n=0; n<nfiles; n++) {
      if (times.find(n)==times.end()) continue;
      for (int i=0; i<nval; i++) {
	for (unsigned j=0; j<nbins; j++) {
	  for (unsigned k=0; k<nbins; k++)
	    out << setw(18) << histo[n][i][j*nbins+k];
	  out << endl;
	}
      }
    }
  }

  MPI_Finalize();

  return 0;
}

