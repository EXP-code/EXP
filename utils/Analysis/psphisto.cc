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
#include <memory>
#include <string>
#include <cmath>

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
#include <foarray.H>
#include <cxxopts.H>
#include <libvars.H>

enum ProjectionType {Cylindrical=1, Spherical=2};
  
int
main(int argc, char **argv)
{
  double rmin, rmax, zcen, zwid;
  int nbins, pbeg, pend, proj, ibeg, iend, iskip;
  string comp, outfile, infile, runtag, fileType, filePrefix;
  bool rlog, logr;

#ifdef DEBUG
  sleep(20);
#endif  
  
  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  cxxopts::Options
    options ("psphisto",
	     "Compute disk potential, force and density profiles from PSP phase-space output files");

  options.add_options()
    ("h,help", "Print this help message")
    ("F,iletype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("P,prefix", "prefix for phase-space files",
     cxxopts::value<std::string>(filePrefix)->default_value("OUT"))
    ("RMIN", "minimum radius for output",
     cxxopts::value<double>(rmin)->default_value("0.0"))
    ("RMAX", "maximum radius for output",
     cxxopts::value<double>(rmax)->default_value("0.1"))
    ("ZCENTER", "disk midplane",
     cxxopts::value<double>(zcen)->default_value("0.0"))
    ("ZWIDTH", "disk halfwidth",
     cxxopts::value<double>(zwid)->default_value("0.05"))
    ("NBINS", "number of bins",
     cxxopts::value<int>(nbins)->default_value("40"))
    ("IBEG", "first PSP index",
     cxxopts::value<int>(ibeg)->default_value("0"))
    ("IEND", "last PSP index",
     cxxopts::value<int>(iend)->default_value("100"))
    ("ISKIP", "skip PSP interval",
     cxxopts::value<int>(iskip)->default_value("1"))
    ("PBEG", "first particle index",
     cxxopts::value<int>(pbeg)->default_value("0"))
    ("PEND", "last particle index",
     cxxopts::value<int>(pend)->default_value("-1"))
    ("LOG", "use logarithmic scaling for radial axis",
     cxxopts::value<bool>(rlog)->default_value("false"))
    ("PROJ", "projection (1=cyl)",
     cxxopts::value<int>(proj)->default_value("1"))
    ("COMP", "component name",
     cxxopts::value<std::string>(comp)->default_value("disk"))
    ("OUTFILE", "filename prefix",
     cxxopts::value<std::string>(outfile)->default_value("histo"))
    ("INFILE", "phase space file prefix",
     cxxopts::value<std::string>(infile)->default_value("OUT"))
    ("RUNTAG", "EXP run tag",
     cxxopts::value<std::string>(runtag)->default_value("run"))
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

  std::vector<std::string> files;
				// Root looks for existence of files
				// with the given tag
  if (myid==0) {
    for (int i=ibeg; i<=iend; i++) {
      std::ostringstream lab;
      lab << infile << "." << runtag << "."
	  << std::setw(5) << std::right << std::setfill('0') << i;
      std::ifstream in(lab.str());
      if (in) files.push_back(lab.str());
      else break;
      std::cout << "." << i << std::flush;
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
      auto l = std::make_unique<char[]>(sz+1);
      MPI_Bcast(l.get(), sz+1, MPI_CHAR, 0, MPI_COMM_WORLD);
      files.push_back(l.get());
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

      PR::PRptr reader = PR::ParticleReader::createReader
	(fileType, {files[n]}, myid, true);

      times[n] = reader->CurrentTime();

      reader->SelectType(comp);

      int icnt = 0;
      for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {

	if (icnt > pbeg) {
	      
	  vector<double> L(3);
	  L[0] = p->mass*(p->pos[1]*p->vel[2] - p->pos[2]*p->vel[1]);
	  L[1] = p->mass*(p->pos[2]*p->vel[0] - p->pos[0]*p->vel[2]);
	  L[2] = p->mass*(p->pos[0]*p->vel[1] - p->pos[1]*p->vel[0]);
	  
	  if (proj==Cylindrical) {
	    if (p->pos[2] >= zcen-zwid && p->pos[2] <= zcen+zwid) {
	      double R = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	      if (rlog) {
		if (R>0.0) {
		  R = log(R);
		  int indx = static_cast<int>(floor( (R - rmin)/dR ));
		  if (indx >=0 && indx<nbins) {
		    histo[n][indx] += p->mass;
		    for (int k=0; k<3; k++) angmom[n][indx*3+k] += L[k];
		  }
		}
	      } else {
		int indx = static_cast<int>(floor( (R - rmin)/dR ));
		if (indx >=0 && indx<nbins) {
		  histo[n][indx] += p->mass;
		  for (int k=0; k<3; k++) angmom[n][indx*3+k] += L[k];
		}
	      }
	    }
	  }
	  else {
	    // if (PROJ==Spherical) {
	    double R = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	    if (rlog) {
	      if (R>0.0) {
		R = log(R);
		int indx = static_cast<int>(floor( (R - rmin)/dR ));
		if (indx >=0 && indx<nbins) {
		  histo[n][indx] += p->mass;
		  for (int k=0; k<3; k++) angmom[n][indx*3+k] += L[k];
		}
	      }
	    } else {
	      int indx = static_cast<int>(floor( (R - rmin)/dR ));
	      if (indx >=0 && indx<nbins) {
		histo[n][indx] += p->mass;
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

