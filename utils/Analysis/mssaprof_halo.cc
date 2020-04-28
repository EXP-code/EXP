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
#include <tuple>
#include <string>
#include <cmath>

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
#include <interp.h>
#include <SphereSL.H>

#include <localmpi.h>
#include <foarray.H>

#include <VtkGrid.H>

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
static  int    OUTR;
  
static  bool VOLUME;
static  bool SURFACE;
static  bool VSLICE;


void write_output(SphereSL& ortho, int indx, int icnt, double time,
		  std::string& file1, std::string& file2, std::string& file3)
{
  unsigned ncnt = 0;
  int noutV = 7, noutS = 7;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  if (indx>=0) sstr << "." << std::setfill('0') << std::setw(5) << indx;
  sstr << "." << std::setfill('0') << std::setw(5) << icnt;

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
	  
	  ortho.all_eval(r, costh, phi, d, d0, p0, p, fr, ft, fp);
	  
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

	  ortho.all_eval(r, costh, phi, d, d0, p0, p, fr, ft, fp);

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


const unsigned int magic_word = 0x5ecede4;

struct CoefElem
{
  int l, m;
  std::vector<double> cos, sin;
};

typedef std::pair<unsigned, unsigned> LM;
typedef std::vector<std::map<double, std::map<LM, CoefElem>>> CoefData;

std::set<LM> LMset;

std::tuple<int, int, CoefData >
get_coefficients(const std::string& coefs)
{
  int lmax=0, nmax;
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
    int LMsize, numT, ngrps;
    in.read((char *)&LMsize,     sizeof(int));
    in.read((char *)&numT,       sizeof(int));
    in.read((char *)&nmax,       sizeof(int));
    in.read((char *)&ngrps,      sizeof(int));
    std::vector<double> times(numT);
    in.read((char *)&times[0],   sizeof(double)*numT);
      
    // Read LM pairs
    //
    for (int i=0; i<LMsize; i++) {
      int LL, MM;
      in.read((char *)&LL, sizeof(int));
      in.read((char *)&MM, sizeof(int));
      LMset.insert({LL, MM});
      lmax = std::max<int>(lmax, LL);
    }

    // Allocate data base
    //
    ret.resize(ngrps);
    for (int p=0; p<ngrps; p++) {
      for (auto t : times) {
	for (auto lm : LMset) {
	  ret[p][t][lm].l = lm.first;
	  ret[p][t][lm].m = lm.second;
	  ret[p][t][lm].cos.resize(nmax, 0);
	  ret[p][t][lm].sin.resize(nmax, 0);
	}
      }
    }

    for (int p=0; p<ngrps; p++) {
      for (auto t : times) {
	for (auto lm : LMset) {
	  for (int n=0; n<nmax; n++) {
	    in.read((char *)&ret[p][t][lm].cos[n], sizeof(double));
	    in.read((char *)&ret[p][t][lm].sin[n], sizeof(double));
	  }
	}
      }
    }

  } else {
    std::cout << "Could not open coefficient file <" << coefs << "> . . . aborting"
	      << std::endl;
    exit(-3);
  }

#if __GNUC__ > 6
  return {lmax, nmax, ret};
#else
  std::tuple<int, int, CoefData > rdat;
  std::get<0>(rdat) = lmax;
  std::get<1>(rdat) = nmax;
  std::get<2>(rdat) = ret;
  return rdat;
#endif

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

  bool DENS, verbose = false, mask = false, All, PCs = false;
  std::string modelfile;
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
    ("PC,p",
     "make rendering for each PC")
    ("All,a",
     po::value<bool>(&All)->default_value(true),
     "make rendering for all PCs")
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
     po::value<std::string>(&outid)->default_value("mssaprof_halo"),
     "Analysis id name")
    ("coeffile",
     po::value<std::string>(&coeffile)->default_value("coef.file"),
     "coefficient file name from exp_mssa")
    ("modfile",
     po::value<std::string>(&modelfile)->default_value("SLGridSph.model"),
     "SL model filename")
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

  if (vm.count("PC")) PCs = true;

  if (vm.count("noCommand")==0) {
    std::string cmdFile = "mssaprof." + outid + ".cmd_line";
    std::ofstream cmd(cmdFile.c_str());
    if (!cmd) {
      std::cerr << "mssaprof: error opening <" << cmdFile
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
  
  if (not PCs and not All) {
    if (myid==0) std::cout << "All output is off . . . exiting" << std::endl;
    exit(0);
  }


  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  std::vector<double> times;

  // ==================================================
  // Read coefficients
  // ==================================================

  auto data = get_coefficients(coeffile);

  // ==================================================
  // Make SL expansion
  // ==================================================

  int lmax     = std::get<0>(data);
  int nmax     = std::get<1>(data);
  auto & coefs = std::get<2>(data);

  SphericalModelTable halo(modelfile);
  SphereSL::mpi = true;
  SphereSL::NUMR = 4000;
  SphereSL ortho(&halo, lmax, nmax);
  
  if (PCs) {

    for (int indx=0; indx<coefs.size(); indx++) {

      std::vector<std::string> outfiles1, outfiles2, outfiles3;
      std::vector<double> T;
      Matrix expcoef;
      
      int count = 0;
      
      for (auto u : coefs[indx]) {
	
	if (count++ % stride) continue;
	
	expcoef.setsize(0, lmax*(lmax+2), 1, nmax);
	expcoef.zero();
	
	int lindx = 0;
	for (int l=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    LM lm = {l, m};
	    if (LMset.find(lm) != LMset.end()) {
	      for (int n=0; n<nmax; n++) 
		expcoef[lindx][n+1] = u.second[lm].cos[n];
	      lindx++;
	      if (m) {
		for (int n=0; n<nmax; n++)
		  expcoef[lindx][n+1] = u.second[lm].sin[n];
		lindx++;
	      }
	    }
	  }
	}
	
	ortho.install_coefs(expcoef);
	
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
	if (outfiles3.size()) writePVD(prefix.str()+".vslice.pvd",  T, outfiles3);
      }
    
    } // PC loop

  } // Output PCs

  // Sum over all PCs
  if (All) {
    std::vector<std::string> outfiles1, outfiles2, outfiles3;
    std::map<double, Matrix> expcoef;
    
    if (myid==0)
      std::cout << "Number of groups in coefficient array:"
		<< coefs.size() << std::endl;

    for (int indx=0; indx<coefs.size(); indx++) {
    
      int count = 0;
      for (auto u : coefs[indx]) {
	
	if (count++ % stride) continue;

				// Find the time in the amp
	double time = u.first;
	std::map<double, Matrix>::iterator expc = expcoef.find(time);

				// Create the coefficient array
	if (expc == expcoef.end()) {
	  expcoef[time].setsize(0, lmax*(lmax+2), 1, nmax);
	  expcoef[time].zero();
	  expc = expcoef.find(time);
	}

	int lindx = 0;
	for (int l=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    LM lm = {l, m};
	    if (LMset.find(lm) != LMset.end()) {
	      for (int n=0; n<nmax; n++)
		expc->second[lindx][n+1] += u.second[lm].cos[n];
	      lindx++;
	      if (m) {
		for (int n=0; n<nmax; n++)
		  expc->second[lindx][n+1] += u.second[lm].sin[n];
		lindx++;
	      }
	    }
	  }
	}
      }
    }

    int count = 0;
    std::vector<double> T;
    for (auto expc : expcoef) {
      double time = expc.first;
      ortho.install_coefs(expc.second);
	
      std::string file1, file2, file3;
    
      if (myid==0) cout << "Writing TOTAL output for  T=" << time
			<< " . . . " << flush;
      write_output(ortho, -1, count++, time, file1, file2, file3);
      //                   ^
      //                   |
      // Will not write index into file name
      //
      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==0) cout << "done" << endl;
    
      if (myid==0) {
	if (file1.size()) outfiles1.push_back(file1);
	if (file2.size()) outfiles2.push_back(file2);
	if (file3.size()) outfiles3.push_back(file3);
      }
      
      T.push_back(time);
    }
    // Time loop

    // Create PVD file
    //
    if (myid==0) {
      std::ostringstream prefix;
      prefix << runtag << ".total";
      if (outfiles1.size()) writePVD(prefix.str()+".volume.pvd",  T, outfiles1);
      if (outfiles2.size()) writePVD(prefix.str()+".surface.pvd", T, outfiles2);
      if (outfiles3.size()) writePVD(prefix.str()+".vslice.pvd",  T, outfiles3);
    }
    
  } // Sum over all PC

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

