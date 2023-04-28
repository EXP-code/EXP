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

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <SphSL.H>

#include <localmpi.H>
#include <cxxopts.H>
#include <foarray.H>
#include <libvars.H>
  
int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  double RMAX, ZMIN, ZMAX;
  int RBINS, IBEG, IEND, ISKIP, PBEG, PEND, ZBINS;
  std::string OUTFILE, INFILE, RUNTAG, CNAME, fileType, filePrefix;
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
    ("CNAME", "Name for gas component",
     cxxopts::value<string>(CNAME)->default_value("gas"))
    ("RUNTAG", "file containing desired indices for PSP output",
     cxxopts::value<string>(RUNTAG)->default_value("run"))
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
  // Do round robin grid assignment of nodes
  // ==================================================

  std::ofstream indx;
  std::ifstream in;

  std::vector<string> files;
				// Root nodes looks for existence of files
  if (myid==0) {
    for (int i=IBEG; i<=IEND; i++) {
      std::ostringstream lab;
      lab << INFILE << "." 
	  << RUNTAG << "." 
	  << setw(5) << right << setfill('0') << i;
      std::ifstream in(lab.str().c_str());
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

      PR::PRptr reader = PR::ParticleReader::createReader
	(fileType, {files[n]}, myid, true);

      times[n] = reader->CurrentTime();
      histo[n] = vector< vector<double> >(nval);
      for (int k=0; k<nval; k++) 
	histo[n][k] = vector<double>(RBINS*ZBINS, 0.0);
      
      
      vector<Particle> particles;

      reader->SelectType(CNAME);

      auto p = reader->firstParticle();
	
      unsigned icnt = 0;

      while (p) {

	if (icnt > PBEG) {
	  double R = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1] );
	  if (p->pos[2] >= ZMIN && p->pos[2] < ZMAX && R < RMAX) {
	    int indR = static_cast<int>(floor( R/dR ));
	    int indZ = static_cast<int>(floor( (p->pos[2] - ZMIN)/dZ ));
	    if (indR >=0 && indR<RBINS &&
		indZ >=0 && indZ<ZBINS ) {
	      histo[n][0][indZ*RBINS+indR] += p->mass;
	      histo[n][1][indZ*RBINS+indR] += p->mass * p->dattrib[0];
	      histo[n][2][indZ*RBINS+indR] += p->mass * p->dattrib[1];
	      histo[n][3][indZ*RBINS+indR] += p->mass * p->dattrib[0] * p->dattrib[1];
	    }
	  }
	}
	
	if (PEND>0 && icnt>PEND) break;
	p = reader->nextParticle();
	icnt++;
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

