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
				// Boost stuff

#include <boost/random/mersenne_twister.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include "Particle.h"
#include <PSP2.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>
#include <VtkGrid.H>
#include <localmpi.H>
#include <foarray.H>

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
boost::mt19937 random_gen;

string OUTFILE;
double RMIN, RMAX, TIME;
int OUTR, NICE, LMAX, NMAX, MMAX, PARTFLAG, L1, L2;
bool ALL, VOLUME, SURFACE, PROBE;

typedef std::shared_ptr<PSP> PSPptr;

void add_particles(PSPptr psp, int& nbods, vector<Particle>& p)
{
  if (myid==0) {

    int nbody = nbods/numprocs;
    int nbody0 = nbods - nbody*(numprocs-1);

				// Send number of bodies to be received
				// by eacn non-root node
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Particle> t(nbody);
    vector<double> val(nbody);

    SParticle *part = psp->GetParticle();
    Particle bod;
    
    //
    // Root's particles
    //
    for (int i=0; i<nbody0; i++) {
      if (part==0) {
	cerr << "Error reading particle [n=" << 0 << ", i=" << i << "]" << endl;
	exit(-1);
      }
      bod.mass = part->mass();
      for (int k=0; k<3; k++) bod.pos[k] = part->pos(k);
      for (int k=0; k<3; k++) bod.vel[k] = part->vel(k);
      p.push_back(bod);

      part = psp->NextParticle();
    }



    //
    // Send the rest of the particles to the other nodes
    //
    for (int n=1; n<numprocs; n++) {
      
      for (int i=0; i<nbody; i++) {
	if (part==0) {
	  cerr << "Error reading particle [n=" 
	       << n << ", i=" << i << "]" << endl;
	  exit(-1);
	}
	t[i].mass = part->mass();
	for (int k=0; k<3; k++) t[i].pos[k] = part->pos(k);
	for (int k=0; k<3; k++) t[i].vel[k] = part->vel(k);
	part = psp->NextParticle();
      }
  
      for (int i=0; i<nbody; i++) val[i] = t[i].mass;
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 11, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 12, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 13, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 14, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 15, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 16, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, n, 17, MPI_COMM_WORLD);
    }

  } else {

    int nbody;
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<Particle> t(nbody);
    vector<double> val(nbody);
				// Get and pack

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 11, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].mass = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 12, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[0] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 13, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[1] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 14, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].pos[2] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 15, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[0] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 16, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[1] = val[i];

    MPI_Recv(&val[0], nbody, MPI_DOUBLE, 0, 17, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].vel[2] = val[i];

    p.insert(p.end(), t.begin(), t.end());
  }

}

void partition(PSPptr psp, std::string& name, vector<Particle>& p)
{
  p.erase(p.begin(), p.end());
  
  PSPstanza* stanza;

  int iok = 1, nbods = 0;

  if (myid==0) {
    stanza = psp->GetNamed(name);
    if (stanza==0) iok = 0;
    nbods = stanza->comp.nbod;
  }

  MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (iok==0) {
    if (myid==0)
      std::cerr << "Could not find component named <" << name << ">"
		<< std::endl;
    MPI_Finalize();
    exit(-1);
  }

  add_particles(psp, nbods, p);
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

      VtkGrid vtk(OUTR, OUTR, OUTR, -RMAX, RMAX, -RMAX, RMAX, -RMAX, RMAX);

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
  std::string MODFILE, runtag, cname, dir("./");
  bool ALL;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  std::ostringstream sout;
  sout << std::string(60, '-') << std::endl
       << "Compute disk potential, force and density profiles from" << std::endl
       << "PSP phase-space output files" << std::endl
       << std::string(60, '-') << std::endl << std::endl
       << "Allowed options";
  
  po::options_description desc(sout.str());
  desc.add_options()
    ("help,h",                                                                          "Print this help message")
    ("OUT",
     "assume original, single binary PSP files as input")
    ("SPL",
     "assume new split binary PSP files as input")
    ("NICE",                po::value<int>(&NICE)->default_value(0),
     "system priority")
    ("RMIN",                po::value<double>(&RMIN)->default_value(0.0),
     "minimum radius for output")
    ("RMAX",                po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("TIME",                po::value<double>(&TIME)->default_value(0.0),
     "Desired time slice")
    ("LMAX",                po::value<int>(&LMAX)->default_value(4),
     "Maximum harmonic order for spherical expansion")
    ("NMAX",                po::value<int>(&NMAX)->default_value(12),
     "Maximum radial order for spherical expansion")
    ("MMAX",                po::value<int>(&MMAX)->default_value(4),
     "Maximum harmonic order")
    ("L1",                  po::value<int>(&L1)->default_value(0),
     "minimum l harmonic")
    ("L2",                  po::value<int>(&L2)->default_value(100),
     "maximum l harmonic")
    ("OUTR",                po::value<int>(&OUTR)->default_value(40),
     "Number of radial points for output")
    ("PROBE",               po::value<bool>(&PROBE)->default_value(true),
     "Make traces along axes")
    ("SURFACE",             po::value<bool>(&SURFACE)->default_value(true),
     "Make equitorial and vertical slices")
    ("VOLUME",              po::value<bool>(&VOLUME)->default_value(false),
     "Make volume for VTK")
    ("ALL",                 po::value<bool>(&ALL)->default_value(false),
     "Compute output for every time slice")
    ("OUTFILE",             po::value<string>(&OUTFILE)->default_value("sphprof"),
     "Filename prefix")
    ("runtag",              po::value<string>(&runtag)->default_value("run0"),
     "Run tag id")
    ("dir,d",               po::value<std::string>(&dir),
     "directory for SPL files")
    ("MODFILE",             po::value<string>(&MODFILE)->default_value("SLGridSph.model"),
     "Halo model file")
    ("beg",                 po::value<int>(&ibeg)->default_value(0),
     "initial frame in sequence")
    ("end",                 po::value<int>(&iend)->default_value(10000),
     "final frame in sequence")
    ("COMP",                po::value<std::string>(&cname)->default_value("stars"),
     "Compute wake for this component name")
    ;
  
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // ==================================================
  // Print help message and exit
  // ==================================================

  if (vm.count("help")) {
    std::cout << std::endl << desc << std::endl;
    return 0;
  }

  bool SPL = false;
  if (vm.count("SPL")) SPL = true;
  if (vm.count("OUT")) SPL = false;

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

  SphericalModelTable halo(MODFILE);
  SphereSL::mpi = true;
  SphereSL::NUMR = 4000;
  SphereSL ortho(&halo, LMAX, NMAX);

  // ==================================================
  // Begin sequence
  // ==================================================

  for (int n=ibeg; n<=iend; n++) {

    int iok = 1;
    std::ostringstream s0;

    if (myid==0) {
      std::ostringstream s0;
      if (SPL) s0 << "SPL.";
      else     s0 << "OUT.";
      s0 << runtag << "." << std::setw(5) << std::setfill('0') << n;

      std::string file = dir + s0.str();
      std::ifstream in(file);
      if (!in) {
	cerr << "Error opening <" << file << ">" << endl;
	iok = 0;
      }
    }
    
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (iok==0) break;

  
    // ==================================================
    // Open frame list
    // ==================================================
    PSPptr psp;

    if (myid==0) {

      if (SPL) psp = PSPptr(new PSPspl(s0.str(), dir, true));
      else     psp = PSPptr(new PSPout(s0.str(), true));

      tnow = psp->CurrentTime();
      cout << "Beginning halo partition [time=" << tnow
	   << ", index=" << n << "] . . . "  << flush;
    }

    vector<Particle> particles;

    partition(psp, cname, particles);
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

