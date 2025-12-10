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
#include <memory>
#include <tuple>
#include <string>
#include <cmath>
#include <set>

                                // System libs
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include "numerical.H"
#include "Particle.h"
#include "interp.H"
#include "SphSL.H"

#include "writePVD.H"
#include "localmpi.H"
#include "foarray.H"

#include "DataGrid.H"
#include "EXPini.H"
#include "libvars.H"
using namespace __EXP__;	// Reference exputils globals

				// Program info string

const std::string overview = "\nCompute disk potential, force and density profiles from\nMSSA reconstructed coefficient files\n";
  
				// Globals
static  string outid;
static  double RMAX;
static  int    OUTR;
  
static  bool VOLUME;
static  bool SURFACE;
static  bool VSLICE;

void write_output(SphSL& ortho, int indx, int icnt, double time,
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

const unsigned int magic_word = 0x5ecede4;

struct CoefElem
{
  int l, m;
  std::vector<double> cos, sin;
};

typedef std::pair<unsigned, unsigned> LM;
typedef std::vector<std::map<double, std::map<LM, CoefElem>>> CoefData;

std::set<LM> LMset;

std::tuple<int, int, CoefData>
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
  std::tuple<int, int, CoefData> rdat;
  std::get<0>(rdat) = lmax;
  std::get<1>(rdat) = nmax;
  std::get<2>(rdat) = ret;
  return rdat;
#endif

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

  bool DENS, verbose = false, mask = false, All, PCs = false;
  std::string modelfile, coeffile, config;
  int stride;

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v,verbose", "verbose output")
    ("b,mask", "blank empty cells")
    ("X,noCommand", "do not save command line")
    ("p,PC", "make rendering for each group or PC")
    ("a,All", "make rendering for all PCs",
     cxxopts::value<bool>(All)->default_value("true"))
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
     cxxopts::value<std::string>(outid)->default_value("mssaprof_halo"))
    ("coeffile", "coefficient file name from exp_mssa",
     cxxopts::value<std::string>(coeffile)->default_value("coef.file"))
    ("modfile", "SL model filename",
     cxxopts::value<std::string>(modelfile)->default_value("SLGridSph.model"))
    ("r,runtag", "runtag for phase space files",
     cxxopts::value<std::string>(runtag)->default_value("run1"))
    ("s,stride", "stride for time output",
     cxxopts::value<int>(stride)->default_value("1"))
    ("T,template", "Write template options file with current and all default values",
     cxxopts::value<string>(config))
    ("f,input", "Input parameter config file",
     cxxopts::value<string>(config))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Write YAML template config file and exit
  //
  if (vm.count("template")) {
    // Write template file
    //
    if (myid==0) SaveConfig(vm, options, "template.yaml");
    MPI_Finalize();
    return 0;
  }

  if (vm.count("help")) {
    if (myid==0) std::cout << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }
 
  // Read parameters fron the YAML config file
  //
  if (vm.count("input")) {
    try {
      vm = LoadConfig(options, config);
    } catch (cxxopts::OptionException& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }
  
  if (vm.count("verbose")) verbose = true;

  if (vm.count("mask"))    mask = true;

  if (vm.count("PC"))      PCs = true;

  if (vm.count("noCommand")==0 and myid==0) {
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

  auto halo = std::make_shared<SphericalModelTable>(modelfile);
  SphSL::mpi = true;
  SphSL::NUMR = 4000;
  SphSL ortho(halo, lmax, nmax);
  
  if (PCs) {

    for (int indx=0; indx<coefs.size(); indx++) {

      std::vector<std::string> outfiles1, outfiles2, outfiles3;
      std::vector<double> T;
      Eigen::MatrixXd expcoef;
      
      int count = 0;
      
      for (auto u : coefs[indx]) {
	
	if (count++ % stride) continue;
	
	expcoef.resize((lmax+1)*(lmax+1), nmax);
	expcoef.setZero();
	
	int lindx = 0;
	for (int l=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    LM lm = {l, m};
	    if (LMset.find(lm) != LMset.end()) {
	      for (int n=0; n<nmax; n++) 
		expcoef(lindx, n) = u.second[lm].cos[n];
	      if (m) {
		for (int n=0; n<nmax; n++)
		  expcoef(lindx+1, n) = u.second[lm].sin[n];
	      }
	    }
	    if (m) lindx += 2;
	    else   lindx += 1;
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
    std::map<double, Eigen::MatrixXd> expcoef;
    
    if (myid==0)
      std::cout << "Number of groups in coefficient array:"
		<< coefs.size() << std::endl;

    for (int indx=0; indx<coefs.size(); indx++) {
    
      int count = 0;
      for (auto u : coefs[indx]) {
	
	if (count++ % stride) continue;

				// Find the time in the amp
	double time = u.first;
	std::map<double, Eigen::MatrixXd>::iterator expc = expcoef.find(time);

				// Create the coefficient array
	if (expc == expcoef.end()) {
	  expcoef[time].resize((lmax+1)*(lmax+1), nmax);
	  expcoef[time].setZero();
	  expc = expcoef.find(time);
	}

	int lindx = 0;
	for (int l=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    LM lm = {l, m};
	    if (LMset.find(lm) != LMset.end()) {
	      for (int n=0; n<nmax; n++)
		expc->second(lindx, n) += u.second[lm].cos[n];
	      if (m) {
		for (int n=0; n<nmax; n++)
		  expc->second(lindx+1, n) += u.second[lm].sin[n];
	      }
	    }
	    if (m) lindx += 2;
	    else   lindx += 1;
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

