/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute VTK slices, and compute volume
 *  for VTK rendering from EXP output
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
 *  MDW 05/08/19
 *
 ***************************************************************************/

				// C++/STL headers
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <cmath>

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include "Particle.h"
#include "Coefs.H"
#include <interp.H>
#include <SphereSL.H>

#include <writePVD.H>
#include <localmpi.H>
#include <cxxopts.H>
#include <foarray.H>

#include <DataGrid.H>

const std::string overview = "Compute disk potential, force and density profiles from\nEXP coefficient files\n";

				// Globals
static  string outid;
static  string runtag;
static  double RMAX;
static  int    OUTR;
  
static  bool VOLUME;
static  bool SURFACE;
static  bool VSLICE;


void write_output(SphereSL& ortho, int indx, double time,
		  std::string& file1, std::string& file2, std::string& file3)
{
  unsigned ncnt = 0;
  int noutV = 7, noutS = 7;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << std::setfill('0') << std::setw(5) << indx;

  string suffix[7] = {"p0", "p", "fr", "ft", "fp", "d0", "d"};

  if (VOLUME) {
      
    // ==================================================
    // Write volume density
    // ==================================================
    
    double dR = 2.0*RMAX/(OUTR-1);
    double x, y, z, r, costh, phi;
    double p0, d0, p, fr, ft, fp, d;
    
    size_t blSiz = OUTR*OUTR*OUTR;
    vector<double> indat(noutV*blSiz, 0.0), otdat(noutV*blSiz);
    
    for (int j=0; j<OUTR; j++) {
	  
      x = -RMAX + dR*j;
	  
      for (int l=0; l<OUTR; l++) {
	
	y = -RMAX + dR*l;
	
	for (int k=0; k<OUTR; k++) {
      
	  z = -RMAX + dR*k;

	  r = sqrt(x*x + y*y + z*z);
	  costh = z/r;
	  phi = atan2(y, x);
	  
	  ortho.all_eval(r, costh, phi, d0, d, p0, p, fr, ft, fp);
	  
	  size_t indx = (OUTR*j + l)*OUTR + k;

	  indat[0*blSiz + indx] = p0;
	  indat[1*blSiz + indx] = p;
	  indat[2*blSiz + indx] = fr;
	  indat[3*blSiz + indx] = ft;
	  indat[4*blSiz + indx] = fp;
	  indat[5*blSiz + indx] = d0;
	  indat[6*blSiz + indx] = d;
	}
      }
    }
    
    
    MPI_Reduce(&indat[0], &otdat[0], noutV*OUTR*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      DataGrid vtk(OUTR, OUTR, OUTR, -RMAX, RMAX, -RMAX, RMAX, -RMAX, RMAX);

      std::vector<double> data(OUTR*OUTR*OUTR);

      for (int n=0; n<noutV; n++) {

	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    for (int k=0; k<OUTR; k++) {
	      size_t indx = (j*OUTR + l)*OUTR + k;

	      data[indx] = otdat[n*blSiz + indx];
	    }
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      std::ostringstream sout;
      sout << runtag +"_" + outid + "_volume" + sstr.str();
      vtk.Write(sout.str());
      file1 = sout.str() + ".vtr";
    }

  }
  
  if (SURFACE) {
    
    // ==================================================
    // Write surface profile
    //   --- in plane ---
    // ==================================================
    
    double dR = 2.0*RMAX/(OUTR-1);
    double x, y, z=0.0, r, costh=0.0, phi;
    double p0, d0, d, p, fr, ft, fp;
    
    vector<double> indat(noutS*OUTR*OUTR, 0.0), otdat(noutS*OUTR*OUTR);
    
    for (int j=0; j<OUTR; j++) {
	
      x = -RMAX + dR*j;
      
      for (int l=0; l<OUTR; l++) {
      
	y = -RMAX + dR*l;
      
	if ((ncnt++)%numprocs == myid) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  ortho.all_eval(r, costh, phi, d0, d, p0, p, fr, ft, fp);
	  
	  indat[(0*OUTR+j)*OUTR+l] = p0;
	  indat[(1*OUTR+j)*OUTR+l] = p;
	  indat[(2*OUTR+j)*OUTR+l] = fr;
	  indat[(3*OUTR+j)*OUTR+l] = ft;
	  indat[(4*OUTR+j)*OUTR+l] = fp;
	  indat[(5*OUTR+j)*OUTR+l] = d0;
	  indat[(6*OUTR+j)*OUTR+l] = d;
	}
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], noutS*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      DataGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::vector<double> data(OUTR*OUTR);

      for (int n=0; n<noutS; n++) {
	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    data[j*OUTR + l] = otdat[(n*OUTR+j)*OUTR+l];
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_surface" + sstr.str();
      vtk.Write(sout.str());
      file2 = sout.str() + ".vtr";
    }
  } // END: SURFACE

  if (VSLICE) {
    
    // ==================================================
    // Write surface profile
    //   --- perp to plane ---
    // ==================================================
    
    double dR = 2.0*RMAX/(OUTR-1);
    double x, y=0, z, r, costh, phi;
    double p0, d0, d, p, fr, ft, fp;
    
    std::vector<double> indat(noutV*OUTR*OUTR, 0.0), otdat(noutV*OUTR*OUTR);
    
    for (int j=0; j<OUTR; j++) {
	
      x = -RMAX + dR*j;

      for (int l=0; l<OUTR; l++) {
      
	z = -RMAX + dR*l;
      
	if ((ncnt++)%numprocs == myid) {
	  
	  r     = sqrt(x*x + z*z);
	  costh = z/r;
	  if (x>=0) phi = 0.0;
	  else      phi = M_PI;

	  ortho.all_eval(r, costh, phi, d0, d, p0, p, fr, ft, fp);

	  indat[(0*OUTR+j)*OUTR+l] = p0;
	  indat[(1*OUTR+j)*OUTR+l] = p;
	  indat[(2*OUTR+j)*OUTR+l] = fr;
	  indat[(3*OUTR+j)*OUTR+l] = ft;
	  indat[(4*OUTR+j)*OUTR+l] = fp;
	  indat[(5*OUTR+j)*OUTR+l] = d0;
	  indat[(6*OUTR+j)*OUTR+l] = d;
	}
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], noutV*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      DataGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::vector<double> data(OUTR*OUTR);

      for (int n=0; n<noutV; n++) {
	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    data[j*OUTR + l] = otdat[(n*OUTR+j)*OUTR+l];
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_vslice" + sstr.str();
      vtk.Write(sout.str());
      file3 = sout.str() + ".vtr";
    }
  } // END: VSLICE

}

std::vector<SphCoefsPtr>
spherical_read(const std::string& file, unsigned stride)
{
  std::vector<SphCoefsPtr> ret;
  std::ifstream in(file);

  unsigned counter = 0;

  while (in.good()) {
    SphCoefsPtr c = std::make_shared<SphCoefs>();
    if (not c->read(in, true)) break;
    if (counter++ % stride == 0) ret.push_back(c);
  }

  return ret;
}


int
main(int argc, char **argv)
{
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  

  //--------------------------------------------------
  // Command-line parsing
  //--------------------------------------------------
  std::string cmd_line;
  for (int i=0; i<argc; i++) {
    cmd_line += argv[i];
    cmd_line += " ";
  }

  bool DENS, verbose = false, mask = false;
  std::string modelfile, coeffile;
  int stride, ibeg, iend;

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v,verbose", "verbose output")
    ("b,mask", "blank empty cells")
    ("X,noCommand", "do not save command line")
    ("R,RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("0.1"))
    ("outr", "number of radial points for output",
     cxxopts::value<int>(OUTR)->default_value("40"))
    ("surface", "make surface equitorial slices",
     cxxopts::value<bool>(SURFACE)->default_value("true"))
    ("vslice", "make surface vertical slices",
     cxxopts::value<bool>(VSLICE)->default_value("false"))
    ("volume", "make volume for VTK rendering",
     cxxopts::value<bool>(VOLUME)->default_value("false"))
    ("o,outid", "Analysis id name",
     cxxopts::value<std::string>(outid)->default_value("halocoef"))
    ("c,coeffile", "coefficient file name",
     cxxopts::value<std::string>(coeffile)->default_value("coef.file"))
    ("m,modfile", "SLGrid model filename",
     cxxopts::value<std::string>(modelfile)->default_value("SLGridSph.model"))
    ("r,runtag", "runtag for phase space files",
     cxxopts::value<std::string>(runtag)->default_value("run1"))
    ("s,stride", "stride for time output",
     cxxopts::value<int>(stride)->default_value("1"))
    ("beg", "initial coefficient frame",
     cxxopts::value<int>(ibeg)->default_value("0"))
    ("end", "final coefficient frame",
     cxxopts::value<int>(iend)->default_value(std::to_string(std::numeric_limits<int>::max())))
    ;
  
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0) std::cout << options.help() << std::endl;
    MPI_Finalize();
    return 1;
  }
 
  if (vm.count("verbose")) verbose = true;

  if (vm.count("mask")) mask = true;

  if (vm.count("noCommand")==0) {
    std::string cmdFile = "haloprof." + outid + ".cmd_line";
    std::ofstream cmd(cmdFile.c_str());
    if (!cmd) {
      std::cerr << "haloprof: error opening <" << cmdFile
		<< "> for writing" << std::endl;
    } else {
      cmd << cmd_line << std::endl;
    }
    
    cmd.close();
  }

#ifdef DEBUG
  sleep(20);
#endif  
  
  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  std::vector<double> times;

  // ==================================================
  // Read coefficients
  // ==================================================

  auto data = spherical_read(coeffile, stride);

  if (data.size()==0) {
    std::cerr << argv[0] << ": no data read from coefficient file <"
	      << coeffile << ">?" << std::endl;
    exit(-1);
  }

  // ==================================================
  // Make SL expansion
  // ==================================================

  int lmax     = data[0]->header.Lmax;
  int nmax     = data[0]->header.nmax;

  auto halo = std::make_shared<SphericalModelTable>(modelfile);
  SphereSL::mpi = true;
  SphereSL::NUMR = 4000;

  SphereSL ortho(halo, lmax, nmax);
  
  std::vector<std::string> outfiles1, outfiles2, outfiles3;
  std::vector<double> T;

  for (int indx=ibeg; indx<=std::min<int>(iend, data.size()); indx++) {
    auto d = data[indx];

    Eigen::MatrixXd expcoef;
    
    int LL = d->header.Lmax + 1;
    expcoef.resize(LL*LL, d->header.nmax);
    expcoef.setZero();
	
    int lindx = 0;
    for (int l=0; l<LL; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<d->header.nmax; n++) 
	  expcoef(lindx, n) = d->coefs(lindx, n);
	if (m) {
	  for (int n=0; n<d->header.nmax; n++)
	    expcoef(lindx+1, n) = d->coefs(lindx+1, n);
	}
	if (m) lindx += 2;
	else   lindx += 1;
      }
    }

    ortho.install_coefs(expcoef);
	
    std::string file1, file2, file3;
	
    if (myid==0) cout << "Writing output for T="
		      << d->header.tnow << " . . . " << flush;
    write_output(ortho, indx, d->header.tnow, file1, file2, file3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    if (myid==0) {
      if (file1.size()) outfiles1.push_back(file1);
      if (file2.size()) outfiles2.push_back(file2);
      if (file3.size()) outfiles3.push_back(file3);
    }
    
    T.push_back(d->header.tnow);
    
  } // Time loop
      
  // Create PVD file
  //
  if (myid==0) {
    std::ostringstream prefix;
    prefix << runtag << "." << ibeg << "_" << iend;
    if (outfiles1.size()) writePVD(prefix.str()+".volume.pvd",  T, outfiles1);
    if (outfiles2.size()) writePVD(prefix.str()+".surface.pvd", T, outfiles2);
    if (outfiles3.size()) writePVD(prefix.str()+".vslice.pvd",  T, outfiles3);
  }
  
  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

