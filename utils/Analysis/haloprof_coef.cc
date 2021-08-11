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
#include <sstream>
#include <tuple>
#include <string>
#include <cmath>

				// BOOST stuff
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp> 

namespace po = boost::program_options;
namespace pt = boost::property_tree;


                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include "Particle.h"
#include "Coefs.H"
#include <interp.H>
#include <SphereSL.H>

#include <localmpi.H>
#include <foarray.H>

#include <VtkGrid.H>

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

      VtkGrid vtk(OUTR, OUTR, OUTR, -RMAX, RMAX, -RMAX, RMAX, -RMAX, RMAX);

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
      
      VtkGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

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
  int stride;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("verbose,v",
     "verbose output")
    ("mask,b",
     "blank empty cells")
    ("noCommand,X",
     "do not save command line")
    ("RMAX,R",
     po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("outr",
     po::value<int>(&OUTR)->default_value(40), 
     "number of radial points for output")
    ("surface",
     po::value<bool>(&SURFACE)->default_value(true),
     "make surface equitorial slices")
    ("vslice",
     po::value<bool>(&VSLICE)->default_value(false),
     "make surface vertical slices")
    ("volume",
     po::value<bool>(&VOLUME)->default_value(false),
     "make volume for VTK rendering")
    ("outid,o",
     po::value<std::string>(&outid)->default_value("halocoef"),
     "Analysis id name")
    ("coeffile,c",
     po::value<std::string>(&coeffile)->default_value("coef.file"),
     "coefficient file name")
    ("modfile,m",
     po::value<std::string>(&modelfile)->default_value("SLGridSph.model"),
     "SLGrid model filename")
    ("runtag,r",
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
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
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

  auto halo = boost::make_shared<SphericalModelTable>(modelfile);
  SphereSL::mpi = true;
  SphereSL::NUMR = 4000;

  SphereSL ortho(halo, lmax, nmax);
  
  std::vector<std::string> outfiles1, outfiles2, outfiles3;
  std::vector<double> T;

  int indx = 0;
  for (auto d: data) {

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
    write_output(ortho, indx++, d->header.tnow, file1, file2, file3);
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
    prefix << runtag << "." << indx;
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

