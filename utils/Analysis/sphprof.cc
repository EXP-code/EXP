/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute gnuplot slices, and compute
 *  volume for VTK rendering
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
 *  MDW 11/28/08, 11/21/19
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
#include <memory>

using namespace std;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>
#include <DataGrid.H>
#include <localmpi.H>
#include <foarray.H>
#include <cxxopts.H>

// Globals

string OUTFILE;
double RMIN, RMAX, TIME;
int OUTR, NICE, LMAX, NMAX, MMAX, PARTFLAG, L1, L2;
bool ALL, VOLUME, SURFACE, PROBE;

void add_particles(PR::PRptr reader, std::string& name, std::vector<Particle>& p)
{
  reader->SelectType(name);

  //
  // Root's particles
  //
  for (auto part=reader->firstParticle(); part != 0; part=reader->nextParticle()) {
    p.push_back(*part);
  }
}

typedef struct {
  double  x;
  double  y;
  double  z;
  double  value;
  int valid;
} Node;


void write_output(SphereSL& ortho, int icnt, double time)
{
  unsigned ncnt = 0;
  Node node;
  int nout;
  
  node.valid = 1;

  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << icnt;

  nout = 8;
  string suffix[8] = {"p0", "p", "fr", "ft", "fp", "d0", "d", "dd"};

  if (VOLUME) {
      
    // ==================================================
    // Write volume density
    // ==================================================
    
    double v;
    int valid = 1;
      
    double dR = 2.0*RMAX/(OUTR-1);
    double x, y, z, r, phi, costh;
    double p0, p1, d0, d1, pl, fr, ft, fp;
    
    vector<double> indat(nout*OUTR*OUTR*OUTR, 0.0), otdat(nout*OUTR*OUTR*OUTR);
    
    for (int k=0; k<OUTR; k++) {
      
      z = -RMAX + dR*k;
      
      for (int l=0; l<OUTR; l++) {
	
	y = -RMAX + dR*l;
	
	for (int j=0; j<OUTR; j++) {
	  
	  x = -RMAX + dR*j;
	  
	  r = sqrt(x*x + y*y + z*z) + 1.0e-18;
	  costh = z/r;
	  phi = atan2(y, x);
	  
	  ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp, L1, L2);
	  
	  indat[((0*OUTR + k)*OUTR + l)*OUTR + j] = p0;
	  indat[((1*OUTR + k)*OUTR + l)*OUTR + j] = p1;
	  indat[((2*OUTR + k)*OUTR + l)*OUTR + j] = fr;
	  indat[((3*OUTR + k)*OUTR + l)*OUTR + j] = ft;
	  indat[((4*OUTR + k)*OUTR + l)*OUTR + j] = fp;
	  indat[((5*OUTR + k)*OUTR + l)*OUTR + j] = d0;
	  indat[((6*OUTR + k)*OUTR + l)*OUTR + j] = d1;
	  if (d0>0.0)
	    indat[((7*OUTR + k)*OUTR + l)*OUTR + j] = d1/d0;
	  else
	    indat[((7*OUTR + k)*OUTR + l)*OUTR + j] = 0.0;
	}
      }
    }
    
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      DataGrid vtk(OUTR, OUTR, OUTR, -RMAX, RMAX, -RMAX, RMAX, -RMAX, RMAX);

      std::vector<double> data(OUTR*OUTR*OUTR);

      for (int n=0; n<nout; n++) {
	for (int k=0; k<OUTR; k++) {
	  for (int l=0; l<OUTR; l++) {
	    for (int j=0; j<OUTR; j++) {
	      data[(j*OUTR + l)*OUTR + k] = otdat[((n*OUTR + k)*OUTR + l)*OUTR + j];
	    }
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      std::ostringstream sout;
      sout << OUTFILE + "_volume";
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
    double x, y, z=0.0, r, phi, costh;
    double p0, p1, d0, d1, fr, ft, fp;
    
    vector<double> indat(nout*OUTR*OUTR, 0.0), otdat(nout*OUTR*OUTR);
    
    for (int l=0; l<OUTR; l++) {
      
      y = -RMAX + dR*l;
      
      for (int j=0; j<OUTR; j++) {
	
	if ((ncnt++)%numprocs == myid) {
	  
	  x = -RMAX + dR*j;
	  
	  r = sqrt(x*x + y*y + z*z) + 1.0e-18;
	  costh = z/r;
	  phi = atan2(y, x);

	  ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	  indat[(0*OUTR+l)*OUTR+j] = p0;
	  indat[(1*OUTR+l)*OUTR+j] = p1;
	  indat[(2*OUTR+l)*OUTR+j] = fr;
	  indat[(3*OUTR+l)*OUTR+j] = ft;
	  indat[(4*OUTR+l)*OUTR+j] = fp;
	  indat[(5*OUTR+l)*OUTR+j] = d0;
	  indat[(6*OUTR+l)*OUTR+j] = d1;
	  if (d0>0.0)
	    indat[(7*OUTR+l)*OUTR+j] = d1/d0;
	  else
	    indat[(7*OUTR+l)*OUTR+j] = 0.0;
	}
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      string name = OUTFILE + ".surf";

      ofstream out(name.c_str());

      if (out) {

	// ==================================================
	// Horizontal line
	// ==================================================
	for (int n=0; n<nout+2; n++)
	  if (n==0) out << "#" << setw(17) << setfill('-') << '-';
	  else out << "|" << setw(17) << setfill('-') << '-';
	out << endl << setfill(' ');
	// ==================================================
	// Field names
	// ==================================================
	out << "# " << setw(16) << left << "x"
	    << "| " << setw(16) << left << "y";
	for (int n=0; n<nout; n++)
	  out << "| " << setw(16) << left << suffix[n];
	out << endl;
	// ==================================================
	// Field index
	// ==================================================
	for (int n=0; n<nout+2; n++)
	  if (n==0) out << "# " << setw(16) << n+1;
	  else out << "| " << setw(16) << n+1;
	out << endl;
	// ==================================================
	// Horizontal line
	// ==================================================
	for (int n=0; n<nout+2; n++)
	  if (n==0) out << "#" << setw(17) << setfill('-') << '-';
	  else out << "|" << setw(17) << setfill('-') << '-';
	out << endl << setfill(' ');
	
	// ==================================================
	// Surface data in GNUPLOT format
	// ==================================================
	for (int l=0; l<OUTR; l++) {
	  y = -RMAX + dR*l;
	
	  for (int j=0; j<OUTR; j++) {
	    x = -RMAX + dR*j;
	    
	    out << setw(18) << x << setw(18) << y;
	    for (int n=0; n<nout; n++)
	      out << setw(18) << otdat[(n*OUTR+l)*OUTR+j];
	    out << endl;
	  }

	  out << endl;

	}

      } else {
	cout << "Error opening surface file <" << name << "> for output"
	     << endl;
      }
    }
  }

  if (PROBE) {
    
    // ==================================================
    // Write line profile along three axes
    // ==================================================
    
    double v;
    float f;
    bool use_log;
    double dR;

    if (RMIN>0.0) {
      use_log = true;
      dR = (log(RMAX) - log(RMIN))/(OUTR-1);
    } else {
      use_log = false;
      dR = RMAX/(OUTR-1);
    }
      
    double r, phi, costh;
    double p0, p1, d0, d1, fr, ft, fp;
    int indx;
    
    vector<double> indat(3*nout*OUTR, 0.0), otdat(3*nout*OUTR);
    
    for (int l=0; l<OUTR; l++) {
      
      r = dR*l;
      if (use_log) r = RMIN*exp(r);
      
      if ((ncnt++)%numprocs == myid) {
	  
	indx = 3*nout*l;

	costh = 0.0;
	phi   = 0.0;
	ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	indat[indx + 0] = p0;
	indat[indx + 1] = p1;
	indat[indx + 2] = fr;
	indat[indx + 3] = ft;
	indat[indx + 4] = fp;
	indat[indx + 5] = d0;
	indat[indx + 6] = d1;
	if (d0>0.0)
	  indat[indx + 7] = d1/d0;
	else
	  indat[indx + 7] = 0.0;

	costh = 0.0;
	phi   = 0.5*M_PI;
	ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	indx += nout;
	indat[indx + 0] = p0;
	indat[indx + 1] = p1;
	indat[indx + 2] = fr;
	indat[indx + 3] = ft;
	indat[indx + 4] = fp;
	indat[indx + 5] = d0;
	indat[indx + 6] = d1;
	if (d0>0.0)
	  indat[indx + 7] = d1/d0;
	else
	  indat[indx + 7] = 0.0;

	costh = 1.0;
	phi   = 0.0;
	ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	indx += nout;
	indat[indx + 0] = p0;
	indat[indx + 1] = p1;
	indat[indx + 2] = fr;
	indat[indx + 3] = ft;
	indat[indx + 4] = fp;
	indat[indx + 5] = d0;
	indat[indx + 6] = d1;
	if (d0>0.0)
	  indat[indx + 7] = d1/d0;
	else
	  indat[indx + 7] = 0.0;

      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], 3*nout*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      vector<string> names(nout);
      for (int i=0; i<nout; i++) {
	names[i] = OUTFILE + "." + suffix[i] + ".cut";
	if (ALL) names[i] += sstr.str();
      }

      foarray out(names, true);

      for (int l=0; l<OUTR; l++) {
	
	r = dR*l;
	if (use_log) r = RMIN*exp(r);
      
	indx = 3*nout*l;
	
	for (int n=0; n<nout; n++)
	  out[n] << setw(18) << time << setw(18) << r
		 << setw(18) << otdat[indx + 0*nout + n]
		 << setw(18) << otdat[indx + 1*nout + n]
		 << setw(18) << otdat[indx + 2*nout + n]
		 << endl;
      }
      
      for (int n=0; n<nout; n++) out[n] << endl;
    }
  }
}


int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  int NICE, LMAX, NMAX, ibeg, iend;
  std::string MODFILE, runtag, cname, dir("./"), fileType;
  bool ALL;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Compute disk potential, force and density profiles from" << std::endl
       << "PSP phase-space output files" << std::endl
       << std::string(60, '-') << std::endl;
  
  cxxopts::Options options(argv[0], sout.str());
  
  options.add_options()
    ("h,help", "Print this help message")
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("NICE", "system priority",
     cxxopts::value<int>(NICE)->default_value("0"))
    ("RMIN", "minimum radius for output",
     cxxopts::value<double>(RMIN)->default_value("0.0"))
    ("RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("0.1"))
    ("TIME", "Desired time slice",
     cxxopts::value<double>(TIME)->default_value("0.0"))
    ("LMAX", "Maximum harmonic order for spherical expansion",
     cxxopts::value<int>(LMAX)->default_value("4"))
    ("NMAX", "Maximum radial order for spherical expansion",
     cxxopts::value<int>(NMAX)->default_value("12"))
    ("MMAX", "Maximum harmonic order",
     cxxopts::value<int>(MMAX)->default_value("4"))
    ("L1", "minimum l harmonic",
     cxxopts::value<int>(L1)->default_value("0"))
    ("L2", "maximum l harmonic",
     cxxopts::value<int>(L2)->default_value("100"))
    ("OUTR", "Number of radial points for output",
     cxxopts::value<int>(OUTR)->default_value("40"))
    ("PROBE", "Make traces along axes",
     cxxopts::value<bool>(PROBE)->default_value("true"))
    ("SURFACE", "Make equitorial and vertical slices",
     cxxopts::value<bool>(SURFACE)->default_value("true"))
    ("VOLUME", "Make volume for VTK",
     cxxopts::value<bool>(VOLUME)->default_value("false"))
    ("ALL", "Compute output for every time slice",
     cxxopts::value<bool>(ALL)->default_value("false"))
    ("OUTFILE", "Filename prefix",
     cxxopts::value<string>(OUTFILE)->default_value("sphprof"))
    ("runtag", "Run tag id",
     cxxopts::value<string>(runtag)->default_value("run0"))
    ("d,dir", "directory for SPL files",
     cxxopts::value<std::string>(dir))
    ("MODFILE", "Halo model file",
     cxxopts::value<string>(MODFILE)->default_value("SLGridSph.model"))
    ("beg", "initial frame in sequence",
     cxxopts::value<int>(ibeg)->default_value("0"))
    ("end", "final frame in sequence",
     cxxopts::value<int>(iend)->default_value("10000"))
    ("COMP", "Compute wake for this component name",
     cxxopts::value<std::string>(cname)->default_value("stars"))
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
    std::cout << std::endl << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }

  if (vm.count("noCommand")==0) {
    std::string cmdFile = runtag + "." + OUTFILE + ".cmd_line";
    std::ofstream cmd(cmdFile);
    if (!cmd) {
      std::cerr << "diskprof4: error opening <" << cmdFile
		<< "> for writing" << std::endl;
    } else {
      std::string cmd_line;
      for (int i=0; i<argc; i++) {
	cmd_line += argv[i];
	cmd_line += " ";
      }
      cmd << cmd_line << std::endl;
    }
    
    cmd.close();
  }

  // ==================================================
  // Check directory for trailing '/'
  // ==================================================
  //
  if (dir.back() != '/') dir += '/';


  // ==================================================
  // Nice process
  // ==================================================

  if (NICE>0)
    setpriority(PRIO_PROCESS, 0, NICE);

  // ==================================================
  // Make SL expansion
  // ==================================================

  auto halo = std::make_shared<SphericalModelTable>(MODFILE);
  SphereSL::mpi = true;
  SphereSL::NUMR = 4000;
  SphereSL ortho(halo, LMAX, NMAX);

  // ==================================================
  // Begin sequence
  // ==================================================

  for (int n=ibeg; n<=iend; n++) {

    auto file0 = PR::ParticleReader::fileNameCreator
      (fileType, n, myid, dir, runtag);

    int iok = 1;
    if (myid==0) {
      std::ifstream in(file0);
      if (!in) {
	cerr << "Error opening <" << file0 << ">" << endl;
	iok = 0;
      }
    }
    
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (iok==0) break;

  
    // ==================================================
    // Open frame list
    // ==================================================
    PR::PRptr reader = PR::ParticleReader::createReader
      (fileType, file0, myid, true);
    
    double tnow = reader->CurrentTime();

    if (myid==0) {
      cout << "Beginning halo partition [time=" << reader->CurrentTime()
	   << ", index=" << n << "] . . . "  << flush;
    }

    std::vector<Particle> particles;

    add_particles(reader, cname, particles);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    if (myid==0) cout << "Accumulating for basis . . . " << flush;
    ortho.reset_coefs();
    for (auto &i : particles) {
      ortho.accumulate(i.pos[0], i.pos[1], i.pos[2], i.mass);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
  
    //------------------------------------------------------------ 

    if (myid==0) cout << "Making coefficients . . . " << flush;
    ortho.make_coefs();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    if (myid==0) cout << "Writing output . . . " << flush;
    write_output(ortho, n, tnow);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
  }

  MPI_Finalize();

  return 0;
}

