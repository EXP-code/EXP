/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute VTK slices, and compute volume
 *  for VTK rendering from MSSA output
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
#include <EmpCylSL.h>

#include <localmpi.h>
#include <foarray.H>

#include <VtkGrid.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

typedef boost::shared_ptr<PSPDump> PSPDumpPtr;

const std::string overview = "Compute disk potential, force and density profiles from\nMSSA reconstructed coefficient files\n";

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
string outdir, runtag, coeffile;
double tpos = 0.0;
double tnow = 0.0;
  
				// Globals
static  string outid;
static  double RMAX;
static  double ZMAX;
static  int    OUTR;
static  int    OUTZ;
  
static  bool VOLUME;
static  bool SURFACE;
static  bool VSLICE;


void write_output(EmpCylSL& ortho, int indx, int icnt, double time,
		  std::string& file1, std::string& file2, std::string& file3)
{
  unsigned ncnt = 0;
  int noutV = 7, noutS = 7;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << std::setfill('0') << std::setw(5) << indx
       << "." << std::setfill('0') << std::setw(5) << icnt;

  string suffix[7] = {"p0", "p", "fr", "fz", "fp", "d0", "d"};

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
    vector<double> indat(noutV*blSiz, 0.0), otdat(noutV*blSiz);
    
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
    
    
    MPI_Reduce(&indat[0], &otdat[0], noutV*OUTZ*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      VtkGrid vtk(OUTR, OUTR, OUTZ, -RMAX, RMAX, -RMAX, RMAX, -ZMAX, ZMAX);

      std::vector<double> data(OUTR*OUTR*OUTZ);

      for (int n=0; n<noutV; n++) {

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
      file1 = sout.str() + ".vtr";
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
    
    vector<double> indat(noutS*OUTR*OUTR, 0.0), otdat(noutS*OUTR*OUTR);
    
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
    
    MPI_Reduce(&indat[0], &otdat[0], noutS*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      VtkGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

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
    
    double v;
    float f;
    
    double dR = 2.0*RMAX/(OUTR-1);
    double dZ = 2.0*ZMAX/(OUTZ-1);
    double x, y=0, z, r, phi;
    double p0, d0, p, fr, fz, fp;
    
    std::vector<double> indat(noutV*OUTR*OUTR, 0.0), otdat(noutV*OUTR*OUTR);
    
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
    
    MPI_Reduce(&indat[0], &otdat[0], noutV*OUTR*OUTZ,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      VtkGrid vtk(OUTR, OUTZ, 1, -RMAX, RMAX, -ZMAX, ZMAX, 0, 0);

      std::vector<double> data(OUTR*OUTZ);

      for (int n=0; n<noutV; n++) {
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
      file3 = sout.str() + ".vtr";
    }
  } // END: VSLICE

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


const unsigned int magic_word = 0x5ecede5;

struct CoefElem
{
  int m;
  std::vector<double> cos, sin;
};

typedef std::vector<std::map<double, CoefElem>> CoefData;

CoefData get_coefficients(const std::string& coefs)
{
  CoefData ret;

  // Open coefficient file
  //
  std::ifstream in(coefs);
  if (in.good()) {

    // Check magic number
    //
    unsigned int test_magic;
    in.read((char *)&test_magic, sizeof(unsigned int));
    if (test_magic != magic_word) {
      std::cout << "Unexpected data in coefficient file <" << coefs << "> . . . aborting"
		<< std::endl;
      exit(-2);
    }

    // Read rest of file
    //
    int MM, numT, nmax, npairs;
    in.read((char *)&MM,         sizeof(int));
    in.read((char *)&numT,       sizeof(int));
    in.read((char *)&nmax,       sizeof(int));
    in.read((char *)&npairs,     sizeof(int));
    std::vector<double> times(numT);
    in.read((char *)&times[0],   sizeof(double)*numT);
      
    // Allocate data base
    //
    ret.resize(npairs);
    for (int p=0; p<npairs; p++) {
      for (auto t : times) {
	ret[p][t].m = MM;
	ret[p][t].cos.resize(nmax, 0);
	ret[p][t].sin.resize(nmax, 0);
      }
    }

    for (int p=0; p<npairs; p++) {
      for (auto t : times) {
	for (int n=0; n<nmax; n++) {
	  in.read((char *)&ret[p][t].cos[n], sizeof(double));
	  in.read((char *)&ret[p][t].sin[n], sizeof(double));
	}
      }
    }

  } else {
    std::cout << "Could not open coefficient file <" << coefs << "> . . . aborting"
	      << std::endl;
    exit(-3);
  }

  return ret;
}

int
main(int argc, char **argv)
{
  int lmax=36, stride=1;
  double rcylmin, rcylmax, rscale, vscale;
  bool DENS, verbose = false, mask = false;
  std::string CACHEFILE;

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
    ("surface",
     po::value<bool>(&SURFACE)->default_value(true),
     "make equatorial slices")
    ("vslice",
     po::value<bool>(&VSLICE)->default_value(true),
     "make vertical slices")
    ("volume",
     po::value<bool>(&VOLUME)->default_value(false),
     "make volume for VTK rendering")
    ("outfile",
     po::value<std::string>(&outid)->default_value("diskprof2"),
     "Filename prefix")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("coeffile",
     po::value<std::string>(&coeffile)->default_value("coef.file"),
     "cachefile name")
    ("runtag",
     po::value<std::string>(&runtag)->default_value("run1"),
     "runtag for phase space files")
    ("stride,s",
     po::value<int>(&stride)->default_value(1), 
     "stride for time output")
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
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  int mmax, numx, numy, nmax, norder, tmp;
  bool cmap=false, dens=false;
  double rmin, rmax, ascl, hscl;

  std::ifstream in(CACHEFILE);
  if (in) {
    in.read((char *)&mmax,   sizeof(int));
    in.read((char *)&numx,   sizeof(int));
    in.read((char *)&numy,   sizeof(int));
    in.read((char *)&nmax,   sizeof(int));
    in.read((char *)&norder, sizeof(int));
    in.read((char *)&tmp,    sizeof(int));    if (tmp) dens = true;
    in.read((char *)&tmp,    sizeof(int));    if (tmp) cmap = true;
    in.read((char *)&rmin,   sizeof(double));
    in.read((char *)&rmax,   sizeof(double));
    in.read((char *)&ascl,   sizeof(double));
    in.read((char *)&hscl,   sizeof(double));
  } else {
    cerr << "mssaprof: error opening cache file <"
	 << CACHEFILE << ">" << std::endl;
    return 0;
  }

  EmpCylSL::RMIN        = rmin;
  EmpCylSL::RMAX        = rmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAP        = cmap;
  EmpCylSL::logarithmic = true;
  EmpCylSL::DENS        = dens;
  EmpCylSL::CACHEFILE   = CACHEFILE;

				// Create expansion
				//
  EmpCylSL ortho(nmax, lmax, mmax, norder, ascl, hscl);
    
  std::vector<double> times;

  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==0) {
    std::cout << "Must have a valid cache  . . aborting" << std::endl;
    exit(-4);
  }

  // ==================================================
  // Read coefficients
  // ==================================================

  auto data = get_coefficients(coeffile);

  for (int indx=0; indx<data.size(); indx++) {

    std::vector<std::string> outfiles1, outfiles2, outfiles3;
    std::vector<double> T;
    int count = 0;

    for (auto u : data[indx]) {

      if (count++ % stride) continue;

      ortho.set_coefs(u.second.m, u.second.cos, u.second.sin, true);

      std::string file1, file2, file3;
    
      if (myid==0) cout << "Writing output for indx=" << indx
			<< ", T=" << u.first << " . . . " << flush;
      write_output(ortho, indx, T.size(), u.first, file1, file2, file3);
      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==0) cout << "done" << endl;
    
      if (myid==0) {
	if (file1.size()) outfiles1.push_back(file1);
	if (file2.size()) outfiles2.push_back(file2);
	if (file3.size()) outfiles3.push_back(file3);
      }
      
      T.push_back(u.first);
    } // Time loop
    
    // Create PVD file
    //
    if (myid==0) {
      std::ostringstream prefix;
	prefix << runtag << "." << indx;
	if (outfiles1.size()) writePVD(prefix.str()+".volume.pvd",  T, outfiles1);
	if (outfiles2.size()) writePVD(prefix.str()+".surface.pvd", T, outfiles2);
	if (outfiles3.size()) writePVD(prefix.str()+".height.pvd",  T, outfiles3);
    }
    
  } // PC loop

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

