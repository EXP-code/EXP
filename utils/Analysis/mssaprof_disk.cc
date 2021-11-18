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
#include <string>
#include <cmath>

#include <yaml-cpp/yaml.h>	// YAML support

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <Particle.H>
#include <interp.H>
#include <EmpCylSL.H>

#include <writePVD.H>
#include <localmpi.H>
#include <cxxopts.H>
#include <foarray.H>

#include <DataGrid.H>

const std::string overview = "Compute disk potential, force and density profiles from\nMSSA reconstructed coefficient files\n";

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

      DataGrid vtk(OUTR, OUTR, OUTZ, -RMAX, RMAX, -RMAX, RMAX, -ZMAX, ZMAX);

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
      
      DataGrid vtk(OUTR, OUTZ, 1, -RMAX, RMAX, -ZMAX, ZMAX, 0, 0);

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

const unsigned int magic_word = 0x5ecede5;

struct CoefElem
{
  std::map<int, std::vector<double>> cos, sin;
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
    int numM, numT, nmin, nmax, ngroups;
    in.read((char *)&numT,       sizeof(int));
    in.read((char *)&nmax,       sizeof(int));
    in.read((char *)&numM,       sizeof(int));
    in.read((char *)&ngroups,    sizeof(int));
    std::vector<int> MM(numM);
    in.read((char *)&MM[0],      sizeof(int)*numM);
    std::vector<double> times(numT);
    in.read((char *)&times[0],   sizeof(double)*numT);
      
    // Allocate data base
    //
    ret.resize(ngroups);
    for (int p=0; p<ngroups; p++) {
      for (auto t : times) {
	for (auto M : MM) {
	  ret[p][t].cos[M].resize(nmax, 0);
	  ret[p][t].sin[M].resize(nmax, 0);
	}
      }
    }

    for (int p=0; p<ngroups; p++) {
      for (auto t : times) {
	for (auto M : MM) {
	  for (int n=0; n<nmax; n++) {
	    in.read((char *)&ret[p][t].cos[M][n], sizeof(double));
	    in.read((char *)&ret[p][t].sin[M][n], sizeof(double));
	  }
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

  int lmax=36, stride=1;
  double rcylmin, rcylmax, rscale, vscale;
  bool DENS, verbose = false, mask = false;
  std::string CACHEFILE, coeffile;

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v,verbose", "verbose output")
    ("b,mask", "blank empty cells")
    ("X,noCommand", "do not save command line")
    ("R,RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("0.1"))
    ("Z,ZMAX", "maximum height for output",
     cxxopts::value<double>(ZMAX)->default_value("0.01"))
    ("outr", "number of radial points for output",
     cxxopts::value<int>(OUTR)->default_value("40"))
    ("outz", "number of vertical points for output",
     cxxopts::value<int>(OUTZ)->default_value("40"))
    ("surface", "make equatorial slices",
     cxxopts::value<bool>(SURFACE)->default_value("true"))
    ("vslice", "make vertical slices",
     cxxopts::value<bool>(VSLICE)->default_value("true"))
    ("volume", "make volume for VTK rendering",
     cxxopts::value<bool>(VOLUME)->default_value("false"))
    ("o,outid", "Analysis id name",
     cxxopts::value<std::string>(outid)->default_value("mssaprof"))
    ("cachefile", "cachefile name",
     cxxopts::value<std::string>(CACHEFILE)->default_value(".eof.cache.file"))
    ("coeffile", "cachefile name",
     cxxopts::value<std::string>(coeffile)->default_value("coef.file"))
    ("r,runtag", "runtag for phase space files",
     cxxopts::value<std::string>(runtag)->default_value("run1"))
    ("s,stride", "stride for time output",
     cxxopts::value<int>(stride)->default_value("1"))
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
    return 0;
  }
 
  if (vm.count("verbose")) verbose = true;

  if (vm.count("mask")) mask = true;

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
  
  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  int mmax, numx, numy, nmax, norder, cmapr=1, cmapz=1, tmp;
  bool dens=false;
  double rmin, rmax, ascl, hscl;

  // Open EOF cachefile
  //
  std::ifstream in(CACHEFILE);

  if (not in) {
    std::cerr << "mssaprof: error opening cache file <"
	      << CACHEFILE << ">" << std::endl;
    return 0;
  }

  // Attempt to read magic number
  //
  unsigned int tmagic;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

  // EOF basis magic number
  //
  const unsigned int hmagic = 0xc0a57a1;

  if (tmagic == hmagic) {
    // YAML size
    //
    unsigned ssize;
    in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));
	
    // Make and read char buffer
    //
    auto buf = std::make_unique<char[]>(ssize+1);
    in.read(buf.get(), ssize);
    buf[ssize] = 0;		// Null terminate
    
    YAML::Node node;
    
    try {
      node = YAML::Load(buf.get());
    }
    catch (YAML::Exception& error) {
      if (myid)
	std::cerr << "YAML: error parsing <" << buf.get() << "> "
		  << "in " << __FILE__ << ":" << __LINE__ << std::endl
		  << "YAML error: " << error.what() << std::endl;
      throw error;
    }
    
    // Get parameters
    //
    mmax    = node["mmax"  ].as<int>();
    numx    = node["numx"  ].as<int>();
    numy    = node["numy"  ].as<int>();
    nmax    = node["nmax"  ].as<int>();
    norder  = node["norder"].as<int>();
    dens    = node["dens"  ].as<bool>();
    if (node["cmap"])
      cmapr = node["cmap"  ].as<int>();
    else 
      cmapr = node["cmapr" ].as<int>();
    if (node["cmapz"])
      cmapz = node["cmapz" ].as<int>();
    rmin    = node["rmin"  ].as<double>();
    rmax    = node["rmax"  ].as<double>();
    ascl    = node["ascl"  ].as<double>();
    hscl    = node["hscl"  ].as<double>();

  } else {
				// Rewind file
    in.clear();
    in.seekg(0);

    in.read((char *)&mmax,   sizeof(int));
    in.read((char *)&numx,   sizeof(int));
    in.read((char *)&numy,   sizeof(int));
    in.read((char *)&nmax,   sizeof(int));
    in.read((char *)&norder, sizeof(int));
    in.read((char *)&dens,   sizeof(int));    if (tmp) dens = true;
    in.read((char *)&cmapr,  sizeof(int));
    in.read((char *)&rmin,   sizeof(double));
    in.read((char *)&rmax,   sizeof(double));
    in.read((char *)&ascl,   sizeof(double));
    in.read((char *)&hscl,   sizeof(double));
  }

  EmpCylSL::RMIN        = rmin;
  EmpCylSL::RMAX        = rmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAPR       = cmapr;
  EmpCylSL::CMAPZ       = cmapz;
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

      bool zero = true;
      for (auto v : u.second.cos) {
	int M = v.first;
	ortho.set_coefs(M, u.second.cos[M], u.second.sin[M], zero);
	zero = false;
      }

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

