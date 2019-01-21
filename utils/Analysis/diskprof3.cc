/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute VTK slices, 
 *  and compute volume for VTK rendering
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
 *  MDW 11/28/08, 02/09/18, 07/11/18
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

				// BOOST stuff
#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp> 

namespace po = boost::program_options;
namespace pt = boost::property_tree;


                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <Vector.h>
#include <numerical.h>
#include "Particle.h"
#include <PSP.H>
#include <interp.h>
#include <EmpOrth9thd.h>

#include <localmpi.h>
#include <foarray.H>

#include <VtkGrid.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

typedef boost::shared_ptr<PSPDump> PSPDumpPtr;

const std::string overview = "Compute disk potential, force and density profiles from\nPSP phase-space output files\n";

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
  
				// Globals
static  string outid;
static  double RMAX;
static  double ZMAX;
static  int OUTR;
static  int OUTZ;
  
static  bool AXIHGT;
static  bool VHEIGHT;
static  bool VOLUME;
static  bool SURFACE;
static  bool VSLICE;

				// Find peak density

double get_max_dens(Vector& vv, double dz)
{
  int lo = vv.getlow();
  int hi = vv.gethigh();

  int ipeak=0;
  double pvalue=-1.0e18;

  for (int i=lo; i<=hi; i++) {
    if (vv[i] > pvalue) {
      ipeak = i;
      pvalue = vv[i];
    }
  }

  if (ipeak == lo) ipeak = lo+1;
  if (ipeak == hi) ipeak = hi-1;

				// Solution of 2nd order Lagrange interpolation
				// formula
  double delta;
  double del = vv[ipeak+1] - vv[ipeak-1];
  double ddel = vv[ipeak+1] + vv[ipeak-1] - 2.0*vv[ipeak];

  if (fabs(ddel) < 1.0e-4) 
    delta = 0.0;		// We're at the peak!
  else
    delta = - 0.5*del/ddel;

  return -ZMAX + dz*(ipeak + delta);

}



Vector get_quart(Vector& vv, double dz)
{
  int lo = vv.getlow();
  int hi = vv.gethigh();

  double next, prev=0.0;
  Vector sum(lo, hi), zz(lo, hi);

  for (int i=lo; i<=hi; i++) {

    if (vv[i] > 0.0)
      next = vv[i];
    else
      next = 0.0;
    
    zz [i] = -ZMAX + dz*(i-lo);
    sum[i] = 0.5*(prev + next);	// Trapezoidal rule
    if (i>lo) sum[i] += sum[i-1];
    prev = next;
  }

  double max = sum[hi];
  Vector ret(-1,1);
  
  ret[-1] = odd2(0.25*max, sum, zz);
  ret[ 0] = odd2(0.50*max, sum, zz);
  ret[ 1] = odd2(0.75*max, sum, zz);

  return ret;
}

Vector get_quart_truncated(Vector& vv, double dz)
{
  int lo = vv.getlow();
  int hi = vv.gethigh();

  int ipeak=0;			// First find peak
  double pvalue=-1.0e18;

  for (int i=lo; i<=hi; i++) {
    if (vv[i] > pvalue) {
      ipeak = i;
      pvalue = vv[i];
    }
  }
				// Zero out above and below first zero
				// from peak

  int lo1 = ipeak;
  int hi1 = ipeak;

  for (; lo1>lo; lo1--) {
    if (vv[lo1]<0.0) break;
  }

  for (; hi1<hi; hi1++) {
    if (vv[hi1]<0.0) break;
  }

  double next, prev=0.0;
  Vector sum(lo1, hi1), zz(lo1, hi1);

  for (int i=lo1; i<=hi1; i++) {

    if (vv[i] > 0.0)
      next = vv[i];
    else
      next = 0.0;
    
    zz [i] = -ZMAX + dz*(i-lo);
    sum[i] = 0.5*(prev + next);	// Trapezoidal rule
    if (i>lo1) sum[i] += sum[i-1];
    prev = next;
  }

  double max = sum[hi1];
  Vector ret(-1,1);
  
  ret[-1] = odd2(0.25*max, sum, zz);
  ret[ 0] = odd2(0.50*max, sum, zz);
  ret[ 1] = odd2(0.75*max, sum, zz);

  return ret;
}


std::string write_output(EmpCylSL& ortho, int icnt, double time)
{
  std::string retstr("nofile");
  unsigned ncnt = 0;
  int nout;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << std::setfill('0') << std::setw(5) << icnt;

  nout = 7;
  string suffix[7] = {"p0", "p", "fr", "fz", "fp", "d0", "d"};

  // ==================================================
  // Axisymmetric structure (for GNUPLOT)
  // ==================================================

  if (AXIHGT) {
    
    double dR = RMAX/(OUTR-1);
    double dz = ZMAX/(OUTZ-1);
    double z, r, phi, hh, d0;
    Vector vv(1, OUTZ);
    Vector q;
    
    vector<double> indat(3*OUTR, 0.0), otdat(3*OUTR);
    
    for (int l=0; l<OUTR; l++) {
      
      if ( (ncnt++)%numprocs == myid ) {
	
	r = dR*l;
	
	for (int k=0; k<OUTZ; k++) {
	  z = -ZMAX + dz*k;
	  phi = 0.0;
	  vv[k+1] = ortho.accumulated_dens_eval(r, z, phi, d0);
	}
	
	indat[0*OUTR+l] = get_max_dens(vv, dz);
	// q = get_quart(vv, dz);
	q = get_quart_truncated(vv, dz);
	indat[1*OUTR+l] = q[0];
	indat[2*OUTR+l] = q[1] - q[-1];
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], 3*OUTR, MPI_DOUBLE, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
    
    if (myid==0) {

      std::string OUTF = runtag + "_" + outid + "_profile" + sstr.str();
      std::ofstream out(OUTF.c_str(), ios::out | ios::app);
      if (!out) {
	cerr << "Error opening <" << OUTF << "> for output\n";
	exit(-1);
      }
      out.setf(ios::scientific);
      out.precision(5);

      for (int l=0; l<OUTR; l++) {
	r = dR*l;
	out
	  << setw(15) << time
	  << setw(15) << r
	  << setw(15) << otdat[0*OUTR+l]
	  << setw(15) << otdat[1*OUTR+l]
	  << setw(15) << otdat[2*OUTR+l]
	  << endl;
      }
      out << endl;
    }
  }

  
  // ==================================================
  // Write vertical position of peak
  // ==================================================
  
  if (VHEIGHT) {
    
    double dR = 2.0*RMAX/(OUTR-1);
    double dz = 2.0*ZMAX/(OUTZ-1);
    double x, y, z, r, phi, hh, d0;
    Vector vv(1, OUTZ);
    Vector q;
    std::vector<double> indat(3*OUTR*OUTR, 0.0), otdat(3*OUTR*OUTR);
    
    for (int j=0; j<OUTR; j++) {
	
      x = -RMAX + dR*j;
	
      for (int l=0; l<OUTR; l++) {
      
	y = -RMAX + dR*l;
      
	if ( (ncnt++)%numprocs == myid ) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  for (int k=0; k<OUTZ; k++)
	    vv[k+1] = ortho.accumulated_dens_eval(r, -ZMAX + dz*k, phi, d0);
	  
	  // q = get_quart(vv, dz);
	  q = get_quart_truncated(vv, dz);
	  
	  indat[(0*OUTR + j)*OUTR + l] =  get_max_dens(vv, dz);
	  indat[(1*OUTR + j)*OUTR + l] =  q[1] - q[-1];
	  indat[(2*OUTR + j)*OUTR + l] =  q[0];
	}
      }
    }
      
      
    MPI_Reduce(&indat[0], &otdat[0], 3*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
      
    if (myid==0) {

      VtkGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::string names[3] = {"surf", "height", "mid"};
      
      std::vector<double> data(OUTR*OUTR);

      for (int i=0; i<3; i++) {

	for (int l=0; l<OUTR; l++) {
	  for (int j=0; j<OUTR; j++) {
	    data[i*OUTR + j] = otdat[(i*OUTR + l)*OUTR + j];
	  }
	}

	vtk.Add(data, names[i]);
      }

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_posn" + sstr.str();
      vtk.Write(sout.str());
    }
  }


  if (VOLUME) {
      
    // ==================================================
    // Write volume density
    // ==================================================
    
    double v;
    int valid = 1;
      
    double dR = 2.0*RMAX/(OUTR-1);
    double dz = 2.0*ZMAX/(OUTZ-1);
    double x, y, z, r, phi;
    double p0, d0, p, fr, fz, fp;
    
    size_t blSiz = OUTZ*OUTR*OUTR;
    vector<double> indat(nout*blSiz, 0.0), otdat(nout*blSiz);
    
    for (int j=0; j<OUTR; j++) {
	  
      x = -RMAX + dR*j;
	  
      for (int l=0; l<OUTR; l++) {
	
	y = -RMAX + dR*l;
	
	for (int k=0; k<OUTZ; k++) {
      
	  z = -ZMAX + dz*k;

	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	  v = ortho.accumulated_dens_eval(r, z, phi, d0);
	  
	  size_t indx = (OUTR*j + l)*OUTR + k;

	  indat[0*blSiz + indx] = p0;
	  indat[1*blSiz + indx] = p;
	  indat[2*blSiz + indx] = fr;
	  indat[3*blSiz + indx] = fz;
	  indat[4*blSiz + indx] = fp;
	  indat[5*blSiz + indx] = d0;
	  indat[6*blSiz + indx] = v;
	}
      }
    }
    
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTZ*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      VtkGrid vtk(OUTR, OUTR, OUTZ, -RMAX, RMAX, -RMAX, RMAX, -ZMAX, ZMAX);

      std::vector<double> data(OUTR*OUTR*OUTZ);

      for (int n=0; n<nout; n++) {

	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    for (int k=0; k<OUTZ; k++) {
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
    }

  }
  
  if (SURFACE) {
    
    // ==================================================
    // Write surface profile
    //   --- in plane ---
    // ==================================================
    
    double v;
    float f;
    
    double dR = 2.0*RMAX/(OUTR-1);
    double x, y, z=0.0, r, phi;
    double p0, d0, p, fr, fz, fp;
    
    vector<double> indat(nout*OUTR*OUTR, 0.0), otdat(nout*OUTR*OUTR);
    
    for (int j=0; j<OUTR; j++) {
	
      x = -RMAX + dR*j;
      
      for (int l=0; l<OUTR; l++) {
      
	y = -RMAX + dR*l;
      
	if ((ncnt++)%numprocs == myid) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	  v = ortho.accumulated_dens_eval(r, z, phi, d0);
	  
	  indat[(0*OUTR+j)*OUTR+l] = p0;
	  indat[(1*OUTR+j)*OUTR+l] = p;
	  indat[(2*OUTR+j)*OUTR+l] = fr;
	  indat[(3*OUTR+j)*OUTR+l] = fz;
	  indat[(4*OUTR+j)*OUTR+l] = fp;
	  indat[(5*OUTR+j)*OUTR+l] = d0;
	  indat[(6*OUTR+j)*OUTR+l] = v;
	}
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      VtkGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::vector<double> data(OUTR*OUTR);

      for (int n=0; n<nout; n++) {
	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    data[l*OUTR + j] = otdat[(n*OUTR+j)*OUTR+l];
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_surface" + sstr.str();
      vtk.Write(sout.str());

      retstr = sout.str();
    }
  } // END: SURFACE

  if (VSLICE) {
    
    // ==================================================
    // Write surface profile
    //   --- perp to plane ---
    // ==================================================
    
    double v;
    float f;
    
    double dR = 2.0*RMAX/(OUTR-1);
    double dZ = 2.0*ZMAX/(OUTZ-1);
    double x, y=0, z, r, phi;
    double p0, d0, p, fr, fz, fp;
    
    std::vector<double> indat(nout*OUTR*OUTR, 0.0), otdat(nout*OUTR*OUTR);
    
      for (int j=0; j<OUTR; j++) {
	
	x = -RMAX + dR*j;

	for (int l=0; l<OUTZ; l++) {
      
	  z = -ZMAX + dZ*l;
      
	  if ((ncnt++)%numprocs == myid) {
	  
	  r   = sqrt(x*x + y*y);
	  phi = atan2(y, x);

	  ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	  v = ortho.accumulated_dens_eval(r, z, phi, d0);
	  
	  indat[(0*OUTR+j)*OUTZ+l] = p0;
	  indat[(1*OUTR+j)*OUTZ+l] = p;
	  indat[(2*OUTR+j)*OUTZ+l] = fr;
	  indat[(3*OUTR+j)*OUTZ+l] = fz;
	  indat[(4*OUTR+j)*OUTZ+l] = fp;
	  indat[(5*OUTR+j)*OUTZ+l] = d0;
	  indat[(6*OUTR+j)*OUTZ+l] = v;
	}
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTZ,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      VtkGrid vtk(OUTR, OUTZ, 1, -RMAX, RMAX, -ZMAX, ZMAX, 0, 0);

      std::vector<double> data(OUTR*OUTZ);

      for (int n=0; n<nout; n++) {
	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTZ; l++) {
	    data[j*OUTZ + l] = otdat[(n*OUTR+j)*OUTZ+l];
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_vslice" + sstr.str();
      vtk.Write(sout.str());
    }
  } // END: VSLICE

  return retstr;
}

void writePVD(const std::string& filename,
	      const std::vector<double>& times,
	      const std::vector<std::string>& files)
{
  // Sanity check
  //
  if (times.size() != files.size()) {
    std::cerr << "Mismatch in file and time arrays" << std::endl;
    exit(-3);
  }

  // Make file collection elements
  //
  pt::ptree ptC;

  for (size_t i=0; i<times.size(); i++) {
    boost::property_tree::ptree x;
    x.put("<xmlattr>.timestep", times[i]);
    x.put("<xmlattr>.part", 0);
    x.put("<xmlattr>.file", files[i]);

    ptC.add_child("DataSet", x);
  }

  // Add VTKFile attributes
  //
  pt::ptree ptP;
  
  ptP.put("<xmlattr>.type", "Collection");
  ptP.put("<xmlattr>.version", "0.1");
  ptP.put("<xmlattr>.byte_order", "LittleEndian");
  ptP.put("<xmlattr>.compressor", "vtkZLibDataCompressor");
  ptP.add_child("Collection", ptC);
  
  // Make the top-level property tree
  //
  pt::ptree PT;

  PT.add_child("VTKFile", ptP);

  // Write the property tree to the XML file.
  //
  pt::xml_parser::write_xml(filename.c_str(), PT, std::locale(), pt::xml_writer_make_settings<std::string>(' ', 4));

  std::cout << "Wrote PVD file <" << filename.c_str() << "> "
	    << " with " << times.size() << " data sets." << std::endl;
}


int
main(int argc, char **argv)
{
  int nice, cnt, nmax, lmax, mmax, nord;
  bool HEIGHT, PVD, verbose = false, mask = false;
  double time, rscale, vscale;
  std::string CACHEFILE, COEFFILE;

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("mask,b",
     "blank empty cells")
    ("nice",
     po::value<int>(&nice)->default_value(0), 
     "number of bins in x direction")
    ("nmax,n",
     po::value<int>(&nmax)->default_value(10), 
     "maximum radial order")
    ("lmax,l",
     po::value<int>(&lmax)->default_value(36), 
     "maximum harmonic order")
    ("mmax,m",
     po::value<int>(&mmax)->default_value(4), 
     "maximum azimuthal order")
    ("norder",
     po::value<int>(&nord)->default_value(4), 
     "maximum empirical radial order")
    ("RMAX,R",
     po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("ZMAX,Z",
     po::value<double>(&ZMAX)->default_value(0.01),
     "maximum height for output")
    ("outr",
     po::value<int>(&OUTR)->default_value(40), 
     "number of radial points for output")
    ("outz",
     po::value<int>(&OUTZ)->default_value(40), 
     "number of vertical points for output")
    ("rscale",
     po::value<double>(&rscale)->default_value(0.01), 
     "radial scale length for basis expansion")
    ("vscale",
     po::value<double>(&vscale)->default_value(0.001), 
     "vertical scale length for basis expansion")
    ("surface",
     po::value<bool>(&SURFACE)->default_value(true),
     "make equatorial slices")
    ("vslice",
     po::value<bool>(&VSLICE)->default_value(true),
     "make vertical slices")
    ("volume",
     po::value<bool>(&VOLUME)->default_value(false),
     "make volume for VTK rendering")
    ("axihgt",
     po::value<bool>(&AXIHGT)->default_value(false),
     "compute midplane height profiles")
    ("height",
     po::value<bool>(&HEIGHT)->default_value(false),
     "compute height profiles")
    ("outfile",
     po::value<std::string>(&outid)->default_value("diskprof2"),
     "Filename prefix")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "basis cache file name")
    ("coeffile",
     po::value<std::string>(&COEFFILE)->default_value("coef.file"),
     "coefficient name")
    ("pvd",
     po::value<bool>(&PVD)->default_value(false),
     "Compute PVD file for ParaView")
    ("runtag",
     po::value<std::string>(&runtag)->default_value("run1"),
     "runtag for phase space files")
    ;
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << overview << std::endl;
    std::cout << desc     << std::endl;
    return 1;
  }
 
  if (vm.count("verbose")) verbose = true;

  if (vm.count("mask")) mask = true;


#ifdef DEBUG
  sleep(20);
#endif  
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  // ==================================================
  // Nice process
  // ==================================================

  if (nice>0) setpriority(PRIO_PROCESS, 0, nice);


  // Create expansion
  //
  EmpCylSL ortho(nmax, lmax, mmax, nord, rscale, vscale);
 
  // Read EOF basis from saved file
  //
  ortho.read_eof_file(CACHEFILE);

  // Set coefficients
  //
  std::ifstream in(COEFFILE);

  std::vector<double> times;
  std::vector<std::string> outfiles;
  
  cnt = 0;

  in >> time;

  while (in) {
    times.push_back(time);

    int M, nmax;
    in >> M;
    in >> nmax;

    Vector cos1(0, nmax-1), sin1(0, nmax-1);
    bool first = true;

    for (int j=0; j<nmax; j++) in >> cos1[j];
    if (M) {
      for (int j=0; j<nmax; j++) in >> sin1[j];
    }

    ortho.set_coefs(M, cos1, sin1, cnt==0);

    outfiles.push_back(write_output(ortho, cnt, time) + ".vtr");

    cnt++;

    in >> time;
  }

  // Create PVD file
  //
  if (myid==0 and PVD) {
    writePVD(runtag+".pvd", times, outfiles);
  }

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

