/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in PSP files for a run and compute spherical or cylindrical
 *  particle distribution
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

using namespace std;

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
  
int
main(int argc, char **argv)
{
  double rmin, rmax, zcen, zwid;
  int nbins, pbeg, pend, proj, ibeg, iend, iskip;
  string comp, outfile, infile, runtag;
  bool rlog, logr;

#ifdef DEBUG
  sleep(20);
#endif  
  
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
    ("RMIN",                po::value<double>(&rmin)->default_value(0.0),
     "minimum radius for output")
    ("RMAX",                po::value<double>(&rmax)->default_value(0.1),
     "maximum radius for output")
    ("ZCENTER",             po::value<double>(&zcen)->default_value(0.0),
     "disk midplane")
    ("ZWIDTH",              po::value<double>(&zwid)->default_value(0.05),
     "disk halfwidth")
    ("NBINS",               po::value<int>(&nbins)->default_value(40),
     "number of bins")
    ("IBEG",                po::value<int>(&ibeg)->default_value(0),
     "first PSP index")
    ("IEND",                po::value<int>(&iend)->default_value(100),
     "last PSP index")
    ("ISKIP",               po::value<int>(&iskip)->default_value(1),
     "skip PSP interval")
    ("PBEG",                po::value<int>(&pbeg)->default_value(0),
     "first particle index")
    ("PEND",                po::value<int>(&pend)->default_value(-1),
     "last particle index")
    ("LOG",                 po::value<bool>(&rlog)->default_value(false),
     "use logarithmic scaling for radial axis")
    ("PROJ",                po::value<int>(&proj)->default_value(1),
     "projection (1=cyl)")
    ("COMP",                po::value<string>(&comp)->default_value("disk"),
     "component name")
     ("LOG",                po::value<bool>(&logr)->default_value(false),
     "use logarithmic scaling for radial axis")
    ("OUTFILE",             po::value<string>(&outfile)->default_value("histo"),
     "filename prefix")
    ("INFILE",              po::value<string>(&infile)->default_value("OUT"),
     "phase space file prefix")
    ("RUNTAG",              po::value<string>(&runtag)->default_value("run"),
     "EXP run tag")
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
				// Root looks for existence of files
				// with the given tag
  if (myid==0) {
    for (int i=ibeg; i<=iend; i++) {
      ostringstream lab;
      lab << infile << "." << runtag << "."
	  << setw(5) << right << setfill('0') << i;
      ifstream in(lab.str().c_str());
      if (in) files.push_back(lab.str());
      else break;
      cout << "." << i << flush;
    }
    cout << endl;
  }

				// Root node sends the file names to
				// everybody
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

  if (rmin>0 && rmax > 0 && rlog) {
    rmin = log(rmin);
    rmax = log(rmax);
  }

  double dR = (rmax - rmin)/(nbins-1);

  vector< vector<double> > histo(nfiles), angmom(nfiles);
  for (int n=0; n<nfiles; n++) {
    histo[n]  = vector<double>(nbins,   0.0);
    angmom[n] = vector<double>(nbins*3, 0.0);
  }

  vector<double> times(nfiles);
  vector<double> rvals(nbins);

				// Set the radial bins
  for (int n=0; n<nbins; n++) {
    rvals[n] = rmin + dR*n;
    if (rlog) rvals[n] = exp(rvals[n]);
  }

  // ==================================================
  // Process the files in parallel
  // ==================================================

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
	  if (stanza->name == comp) break;
	}

	if (stanza==0) {
	  std::cout << "Could not find Component <" << comp << "> at time = "
		    << psp->CurrentTime() << std::endl;
	} else {

	  int icnt = 0;
	  for (SParticle* 
		 p=psp->GetParticle(); p!=0; p=psp->NextParticle()) {

	    if (icnt > pbeg) {
	      
	      vector<double> L(3);
	      L[0] = p->mass()*(p->pos(1)*p->vel(2) - p->pos(2)*p->vel(1));
	      L[1] = p->mass()*(p->pos(2)*p->vel(0) - p->pos(0)*p->vel(2));
	      L[2] = p->mass()*(p->pos(0)*p->vel(1) - p->pos(1)*p->vel(0));

	      if (proj==Cylindrical) {
		if (p->pos(2) >= zcen-zwid && p->pos(2) <= zcen+zwid) {
		  double R = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
		  if (rlog) {
		    if (R>0.0) {
		      R = log(R);
		      int indx = static_cast<int>(floor( (R - rmin)/dR ));
		      if (indx >=0 && indx<nbins) {
			histo[n][indx] += p->mass();
			for (int k=0; k<3; k++) angmom[n][indx*3+k] += L[k];
		      }
		    }
		  } else {
		    int indx = static_cast<int>(floor( (R - rmin)/dR ));
		    if (indx >=0 && indx<nbins) {
		      histo[n][indx] += p->mass();
		      for (int k=0; k<3; k++) angmom[n][indx*3+k] += L[k];
		    }
		  }
		}
	      }
	      else {
		// if (PROJ==Spherical) {
		double R = sqrt(p->pos(0)*p->pos(0) + p->pos(1)*p->pos(1));
		if (rlog) {
		  if (R>0.0) {
		    R = log(R);
		    int indx = static_cast<int>(floor( (R - rmin)/dR ));
		    if (indx >=0 && indx<nbins) {
		      histo[n][indx] += p->mass();
		      for (int k=0; k<3; k++) angmom[n][indx*3+k] += L[k];
		    }
		  }
		} else {
		  int indx = static_cast<int>(floor( (R - rmin)/dR ));
		  if (indx >=0 && indx<nbins) {
		    histo[n][indx] += p->mass();
		    for (int k=0; k<3; k++) angmom[n][indx*3+k] += L[k];
		  }
		}
	      }
	    }
	    
	    if (pend>0 && icnt>pend) break;
	    icnt++;
	  }
	}
      }
    }
  }
  
  if (myid==0) {
    MPI_Reduce(MPI_IN_PLACE, &times[0], nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int n=0; n<nfiles; n++) {
      MPI_Reduce(MPI_IN_PLACE, &histo[n][0], nbins, 
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &angmom[n][0], nbins*3, 
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Reduce(&times[0], 0, nfiles,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int n=0; n<nfiles; n++) {
      MPI_Reduce(&histo[n][0], 0, nbins, 
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&angmom[n][0], 0, nbins*3, 
		 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }

  if (myid==0) {
    ofstream out(outfile + ".dat");

    // Label header
    out << left 
	<< setw(18) << "# Time"
	<< setw(18) << "Radius"
	<< setw(18) << "Mass"
	<< setw(18) << "Cumulated mass"
	<< setw(18) << "Lx/mass"
	<< setw(18) << "Ly/mass"
	<< setw(18) << "Lz/mass"
	<< setw(18) << "Ltot/mass"
	<< setw(18) << "Cumulated Lx"
	<< setw(18) << "Cumulated Ly"
	<< setw(18) << "Cumulated Lz"
	<< endl
	<< setw(18) << "# 1"
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
	<< endl << right;

    for (int n=0; n<nfiles; n++) {
      double sum = 0.0, LL;
      vector<double> Lsum(4, 0.0);
      for (unsigned j=0; j<nbins; j++) {
	out << setw(18) << times[n] << setw(18) << rvals[j]
	    << setw(18) << histo[n][j]
	    << setw(18) << (sum += histo[n][j]);
	LL = 0.0;
	for (unsigned k=0; k<3; k++) {
	  if (histo[n][j]>0.0)
	    out << setw(18) << angmom[n][j*3+k]/histo[n][j];
	  else
	    out << setw(18) << 0.0;
	  Lsum[k] += angmom[n][j*3+k];
	  LL += angmom[n][j*3+k]*angmom[n][j*3+k];
	}
	LL = sqrt(LL);
	if (histo[n][j]>0.0)
	  out << setw(18) << LL/histo[n][j];
	else
	  out << setw(18) << 0.0;
	Lsum[3] += sqrt(LL);
	if (sum>0.0) {
	  for (unsigned k=0; k<3; k++)
	    out << setw(18) << Lsum[k]/sum;
	  out << setw(18) << Lsum[3]/sum;
	} else {
	  for (unsigned k=0; k<4; k++)
	    out << setw(18) << 0.0;
	}
	out << endl;
      }
      out << endl;
    }
    
  }

  MPI_Finalize();

  return 0;
}

