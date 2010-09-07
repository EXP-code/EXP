/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute gnuplot slices, 
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
 *  MDW 11/28/08
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

using namespace std;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <Vector.h>
#include <numerical.h>
#include <Particle.h>
#include <PSP.H>
#include <interp.h>
#include <massmodel.h>
#include <SphereSL.H>

#include <localmpi.h>
#include <ProgramParam.H>
#include <foarray.H>

program_option init[] = {
  {"NICE",		"int",		"0",		"system priority"},
  {"DENS",		"bool",		"true",		"compute density"},
  {"RMIN",		"double",	"0.0",		"minimum radius for output"},
  {"RMAX",		"double",	"0.1",		"maximum radius for output"},
  {"TIME",		"double",	"0.0",		"Desired time slice"},
  {"LMAX",		"int",		"4",		"Maximum harmonic order for spherical expansion"},
  {"NMAX",		"int",		"12",		"Maximum radial order for spherical expansion"},
  {"MMAX",		"int",		"4",		"Maximum harmonic order"},
  {"OUTR",		"int",		"40",		"Number of radial points for output"},
  {"PROBE",		"bool",		"true",		"Make traces along axes"},
  {"SURFACE",		"bool",		"true",		"Make equitorial and vertical slices"},
  {"VOLUME",		"bool",		"false",	"Make volume for VTK"},
  {"ALL",		"bool",		"false",	"Compute output for every time slice"},
  {"PARTFLAG",		"int",		"4",		"Wakes using Component(s) [1=stars | 2=gas | 4=halo]"},
  {"OUTFILE",		"string",	"haloprof",	"Filename prefix"},
  {"MODFILE",		"string",	"SLGridSph.model",
   							"Model file"},
  {"INFILE",		"string",	"OUT",		"Phase space file"},
  {"INDEX",		"string",	"frame.indx",	"File containing desired indices for PSP output"},
  {"",			"",		"",		""}
};


const char desc[] = "Compute disk potential, force and density profiles from PSP phase-space output files\n";


ProgramParam config(desc, init);

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
double tpos = 0.0;
double tnow = 0.0;
  
enum ComponentType {Star=1, Gas=2, Halo=4};

void add_particles(ifstream* in, PSPDump* psp, int& nbods, vector<Particle>& p)
{
  if (myid==0) {

    int nbody = nbods/numprocs;
    int nbody0 = nbods - nbody*(numprocs-1);

				// Send number of bodies to be received
				// by eacn non-root node
    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Particle> t(nbody);
    vector<double> val(nbody);

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
      bod.mass = part->mass;
      for (int k=0; k<3; k++) bod.pos[k] = part->pos[k];
      for (int k=0; k<3; k++) bod.vel[k] = part->vel[k];

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
	t[i].mass = part->mass;
	for (int k=0; k<3; k++) t[i].pos[k] = part->pos[k];
	for (int k=0; k<3; k++) t[i].vel[k] = part->vel[k];

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

void partition(ifstream* in, PSPDump* psp, int cflag, vector<Particle>& p)
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

  if (cflag & Halo) {
    if (myid==0) {
      int nbods = psp->CurrentDump()->ndark;
      psp->GetDark();

      add_particles(in, psp, nbods, p);
    } else {
      add_particles(in, psp, nbods, p);
    }      
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

  string OUTFILE  = config.get<string>("OUTFILE");
  double RMIN     = config.get<double>("RMIN");
  double RMAX     = config.get<double>("RMAX");
  int OUTR        = config.get<int   >("OUTR");
  
  bool ALL        = config.get<bool>("ALL");
  bool VOLUME     = config.get<bool>("VOLUME");
  bool SURFACE    = config.get<bool>("SURFACE");
  bool PROBE      = config.get<bool>("PROBE");

  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << icnt;

  nout = 7;
  string suffix[7] = {"p0", "p", "fr", "ft", "fp", "d0", "d"};

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
	  
	  ortho.all_eval(r, costh, phi, d0, d1, p0, p1, fr, ft, fp);
	  
	  indat[((0*OUTR + k)*OUTR + l)*OUTR + j] = p0;
	  indat[((1*OUTR + k)*OUTR + l)*OUTR + j] = p1;
	  indat[((2*OUTR + k)*OUTR + l)*OUTR + j] = fr;
	  indat[((3*OUTR + k)*OUTR + l)*OUTR + j] = ft;
	  indat[((4*OUTR + k)*OUTR + l)*OUTR + j] = fp;
	  indat[((5*OUTR + k)*OUTR + l)*OUTR + j] = d0;
	  indat[((6*OUTR + k)*OUTR + l)*OUTR + j] = d1;
	}
      }
    }
    
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      vector<string> names(nout);
      
      for (int i=0; i<nout; i++) {
	names[i] = OUTFILE + "." + suffix[i] + ".vol";
	if (ALL) names[i] += sstr.str();
      }

      foarray out(names);

      for (int i=0; i<nout; i++) {
	out[i].write((char *)&OUTR, sizeof(int));
	out[i].write((char *)&OUTR, sizeof(int));
	out[i].write((char *)&OUTR, sizeof(int));
      }
    
    
      for (int k=0; k<OUTR; k++) {
	
	node.z = -RMAX + dR*k;
	
	for (int l=0; l<OUTR; l++) {
	  
	  node.y = -RMAX + dR*l;
	  
	  for (int j=0; j<OUTR; j++) {
	    
	    node.x = -RMAX + dR*j;
	    
	    for (int n=0; n<nout; n++) {
	      node.value  = otdat[((n*OUTR + k)*OUTR + l)*OUTR + j];

	      out[n].write((char *)&node, sizeof(Node));
	    }
	  }
	}
      }
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
  
  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  if (config.parse_args(argc,argv)) return -1;
  
  
  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  // ==================================================
  // Nice process
  // ==================================================

  if (config.get<int>("NICE")>0)
    setpriority(PRIO_PROCESS, 0, config.get<int>("NICE"));


  // ==================================================
  // Make SL expansion
  // ==================================================

  SphericalModelTable halo(config.get<string>("MODFILE"));
  SphereSL::mpi = true;
  SphereSL::NUMR = 4000;
  SphereSL ortho(&halo, config.get<int>("LMAX"), config.get<int>("NMAX"));

  // ==================================================
  // Open frame list
  // ==================================================

  ofstream indx;
  ifstream in;

  if (myid==0) {

    indx.open(config.get<string>("INDEX").c_str());
    if (!indx) {
      cerr << "Error opening <" << config.get<string>("INDEX") 
	   << "> for output . . . continuing without writing index file" 
	   << endl;
    }


    in.open(config.get<string>("INFILE").c_str());
    if (!in) {
      cerr << "Error opening <" << config.get<string>("INFILE")
	   << ">" << endl;
      exit(-1);
    }
  }

  
  PSPDump psp(&in, true);

  Dump *dump = psp.GetDump();
  bool ALL = config.get<bool>("ALL");

  if (!ALL) dump = psp.CurrentDump();

  int icnt = 0;
  vector<Particle> particles;

  while (dump) {

    //------------------------------------------------------------ 

    if (myid==0) {
      cout << "Beginning particle partition [time="
	   << dump->header.time << "] . . . " << flush;
      if (ALL) 
	indx << setw(15) << icnt << setw(15) << dump->header.time << endl;
    }

    if (in.rdstate() & ios::eofbit) {
      in.close();
      in.open(config.get<string>("INFILE").c_str());
    }

    partition(&in, &psp, config.get<int>("PARTFLAG"), particles);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    if (myid==0) cout << "Accumulating for basis . . . " << flush;
    ortho.reset_coefs();
    for (vector<Particle>::iterator it=particles.begin(); it!=particles.end(); it++) 
      {
	ortho.accumulate(it->pos[0], it->pos[1], it->pos[2], it->mass);
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
    write_output(ortho, icnt++, dump->header.time);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;

    //------------------------------------------------------------ 

    dump = psp.NextDump();
  }

  MPI_Finalize();

  return 0;
}

