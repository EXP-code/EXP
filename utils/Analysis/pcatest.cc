/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Compute pca basis and create VTK slices
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
 *  MDW 03/10/18
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
#include <ProgramParam.H>
#include <foarray.H>

#include <VtkGrid.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

typedef boost::shared_ptr<PSPDump> PSPDumpPtr;

const std::string overview = "Compute and print PCA basis for rendering";

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

enum ComponentType {Star=1, Gas=2};

void add_particles(ifstream* in, PSPDumpPtr psp, int& nbods, vector<Particle>& p)
{
  if (myid==0) {

    int nbody = nbods/numprocs;
    int nbody0 = nbods - nbody*(numprocs-1);

				// Send number of bodies to be received
				// by eacn non-root node
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Particle>        t(nbody);
    vector<double>        val(nbody);
    vector<unsigned long> seq(nbody);

    SParticle *part = psp->GetParticle(in);
    Particle bod;
    
    //
    // Root's particles
    //
    for (int i=0; i<nbody0; i++) {
      if (part==0) {
	cerr << "Error reading particle [n=" << 0 << ", i=" << i << "]" << endl;
	exit(-1);
      }

      // Make a Particle
      //
      bod.mass = part->mass();
      for (int k=0; k<3; k++) bod.pos[k] = part->pos(k);
      for (int k=0; k<3; k++) bod.vel[k] = part->vel(k);
      bod.indx = part->indx();
      p.push_back(bod);

      part = psp->NextParticle(in);
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
	part = psp->NextParticle(in);
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

      for (int i=0; i<nbody; i++) seq[i] = t[i].indx;
      MPI_Send(&seq[0], nbody, MPI_UNSIGNED_LONG, n, 18, MPI_COMM_WORLD);
    }

  } else {

    int nbody;
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<Particle>        t(nbody);
    vector<double>        val(nbody);
    vector<unsigned long> seq(nbody);
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

    MPI_Recv(&seq[0], nbody, MPI_UNSIGNED_LONG, 0, 18, 
	     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i=0; i<nbody; i++) t[i].indx = seq[i];

    p.insert(p.end(), t.begin(), t.end());
  }
}


void partition(ifstream* in, PSPDumpPtr psp, int cflag, vector<Particle>& p)
{
  p.erase(p.begin(), p.end());

  int nbods = 0;
  if (cflag & Star) {
    if (myid==0) {
      nbods = psp->CurrentDump()->nstar;
      psp->GetStar();

      add_particles(in, psp, nbods, p);
    } else {
      add_particles(in, psp, nbods, p);
    }      
  }

  if (cflag & Gas) {
    if (myid==0) {
      int nbods = psp->CurrentDump()->ngas;
      psp->GetGas();

      add_particles(in, psp, nbods, p);
    } else {
      add_particles(in, psp, nbods, p);
    }      
  }

}

int
main(int argc, char **argv)
{
  int nice, numx, numy, lmax, mmax, nmax, norder;
  int initc, partc, beg, end, stride, init;
  double rcylmin, rcylmax, rscale, vscale;
  bool verbose = false, DENS;
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
    ("nice",
     po::value<int>(&nice)->default_value(0), 
     "number of bins in x direction")
    ("RMAX,R",
     po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("ZMAX,Z",
     po::value<double>(&ZMAX)->default_value(0.01),
     "maximum height for output")
    ("rcylmin",
     po::value<double>(&rcylmin)->default_value(0.001),
     "minimum radius for cylindrical basis table")
    ("rcylmax",
     po::value<double>(&rcylmax)->default_value(20.0),
     "maximum radius for cylindrical basis table")
    ("NUMX",
     po::value<int>(&numx)->default_value(128), 
     "number of radial table entries")
    ("NUMY",
     po::value<int>(&numy)->default_value(64), 
     "number of vertical table entries")
    ("rscale",
     po::value<double>(&rscale)->default_value(0.01), 
     "radial scale length for basis expansion")
    ("vscale",
     po::value<double>(&vscale)->default_value(0.001), 
     "vertical scale length for basis expansion")
    ("lmax",
     po::value<int>(&lmax)->default_value(36), 
     "maximum harmonic order for spherical expansion")
    ("nmax",
     po::value<int>(&nmax)->default_value(8),
     "maximum harmonic order for spherical expansion")
    ("mmax",
     po::value<int>(&mmax)->default_value(4), 
     "maximum azimuthal harmonic order for cylindrical expansion")
    ("norder",
     po::value<int>(&norder)->default_value(4), 
     "maximum radial order for each harmonic subspace")
    ("outr",
     po::value<int>(&OUTR)->default_value(40), 
     "number of radial points for output")
    ("outz",
     po::value<int>(&OUTZ)->default_value(40), 
     "number of vertical points for output")
    ("density",
     po::value<bool>(&DENS)->default_value(true),
     "compute density")
    ("partflag",
     po::value<int>(&partc)->default_value(1),
     "Wakes using Component(s) [1=stars | 2=gas]")
    ("initcomp",
     po::value<int>(&initc)->default_value(1),
     "train on Component (1=stars)")
    ("init",
     po::value<int>(&init)->default_value(0),
     "fiducial PSP index")
    ("beg",
     po::value<int>(&beg)->default_value(0),
     "initial PSP index")
    ("end",
     po::value<int>(&end)->default_value(99999),
     "final PSP index")
    ("stride",
     po::value<int>(&stride)->default_value(1),
     "PSP index stride")
    ("outfile",
     po::value<std::string>(&outid)->default_value("diskprof2"),
     "Filename prefix")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("runtag",
     po::value<std::string>(&runtag)->default_value("run1"),
     "runtag for phase space files")
    ;
  
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("verbose")) verbose = true;

#ifdef DEBUG
  sleep(20);
#endif  
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  if (vm.count("help")) {
    if (myid==0) {
      std::cout << overview << std::endl;
      std::cout << desc     << std::endl;
    }

    MPI_Finalize();
    return 0;
  }
 
  // ==================================================
  // Nice process
  // ==================================================

  if (nice>0) setpriority(PRIO_PROCESS, 0, nice);


  // ==================================================
  // PSP input stream
  // ==================================================

  int iok = 1;
  ifstream in0, in1;
  std::ostringstream s0, s1;
  if (myid==0) {
    s0 << "OUT." << runtag << "."
       << std::setw(5) << std::setfill('0') << init;
    in0.open(s0.str());
    if (!in0) {
      cerr << "Error opening <" << s0.str() << ">" << endl;
      iok = 0;
    }
  }
    
  MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (iok==0) {
    MPI_Finalize();
    exit(-1);
  }
    
  // ==================================================
  // All processes will now compute the basis functions
  // *****Using MPI****
  // ==================================================

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAP        = true;
  EmpCylSL::logarithmic = true;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::CACHEFILE   = CACHEFILE;

				// Create expansion
				//
  EmpCylSL ortho(nmax, lmax, mmax, norder, rscale, vscale);
    
  vector<Particle> particles;
  PSPDumpPtr psp;
  
  std::vector<double> times;
  std::vector<std::string> outfiles;

  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==0) {
    
    //------------------------------------------------------------ 

    if (myid==0) {
      psp = PSPDumpPtr(new PSPDump (&in0, true));
      cout << "Beginning disk partition [time="
	   << psp->CurrentTime()
	   << "] . . . " << flush;
    }
      
    // Do we need to close and reopen?
    if (in0.rdstate() & ios::eofbit) {
      in0.close();
      in0.open(s0.str());
    }
      
    partition(&in0, psp, initc, particles);
    if (myid==0) cout << "done" << endl;

    if (myid==0) {
      cout << endl
	   << setw(4) << "#" << "  " << setw(8) << "Number" << endl 
	   << setfill('-')
	   << setw(4) << "-" << "  " << setw(8) << "-" << endl
	   << setfill(' ');
    }
    for (int n=0; n<numprocs; n++) {
      if (n==myid) {
	cout << setw(4) << myid << "  "
	     << setw(8) << particles.size() << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if (myid==0) cout << endl;

    //------------------------------------------------------------ 

    if (myid==0) cout << "Beginning disk accumulation . . . " << flush;
    ortho.setup_eof();
    ortho.setup_accumulation();
    ortho.accumulate_eof(particles, true);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
    
    if (myid==0) cout << "Making the EOF . . . " << flush;
    ortho.make_eof();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
    
    if (myid==0) cout << "Making disk coefficients . . . " << flush;
    ortho.make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
  }


  
  // ==================================================
  // Set up for Hall
  // ==================================================
  EmpCylSL::SELECT = true;
  ortho.setHall(runtag + "_" + outid + ".pca", 1);


  // ==================================================
  // The PSP loop
  // ==================================================
  for (int indx=beg; indx<=end; indx+=stride) {

    // ==================================================
    // PSP input stream
    // ==================================================

    iok = 1;
    if (myid==0) {
      s1.str("");		// Clear
      s1 << "OUT." << runtag << "."
	 << std::setw(5) << std::setfill('0') << indx;
      
      in1.close();		// Make sure previous file is closed
      in1.open(s1.str());	// Now, try to open a new one . . . 
      if (!in1) {
	cerr << "Error opening <" << s1.str() << ">" << endl;
	iok = 0;
      }
    }
    
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (iok==0) break;
    
    // ==================================================
    // Open frame list
    // ==================================================
    
    psp = PSPDumpPtr(new PSPDump (&in1, true));

    Dump *dump = psp->GetDump();
    
    dump = psp->CurrentDump();
    
    if (in1.rdstate() & ios::eofbit) {
      in1.close();
      in1.open(s1.str());
    }
      
    //------------------------------------------------------------ 

    partition(&in1, psp, partc, particles);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    if (myid==0) cout << "Accumulating particle positions . . . " << flush;
    ortho.accumulate(particles, 0, true);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
      
    //------------------------------------------------------------ 
      
    if (myid==0) cout << "Making disk coefficients . . . " << flush;
    ortho.make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
      
    //------------------------------------------------------------ 
      
    if (myid==0) {
      cout << "Writing output: (m, n)=" << flush;

      for (int mm=0; mm<=mmax; mm++) {
	for (int nn=1; nn<=norder; nn++) {
	  std::cout << "(" << mm << ", " << nn << ") " << flush;
	  ortho.dump_images_basis_pca(runtag,
				      RMAX, ZMAX, OUTR, OUTZ,
				      mm, nn, indx);
	}
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
    
  } // Dump loop

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

