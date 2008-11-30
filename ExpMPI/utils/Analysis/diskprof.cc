// This may look like C code, but it is really -*- C++ -*-

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
#include <gaussQ.h>
#include <Particle.h>
#include <PSP.H>
#include <interp.h>
#include <EmpOrth9thd.h>

#include <localmpi.h>
#include <ProgramParam.H>

program_option init[] = {
  {"NICE",		"int",		"0",		"system priority"},
  {"DENS",		"bool",		"true",		"compute density"},
  {"RMAX",		"double",	"1.0",		"maximum radius for output"},
  {"ZMAX",		"double",	"1.0",		"maximum height for output"},
  {"RCYLMIN",		"double",	"0.001",	"number of scale lengths for minimum radius in table"},
  {"RCYLMAX",		"double",	"20.0",		"number of scale lengths for maximum radius in table"},
  {"NUMX",		"int",		"128",		"number of radial table entries"},
  {"NUMY",		"int",		"64",		"number of vertical table entries"},
  {"EXPSCALE",		"double",	"0.01",		"Radial scale length for basis expansion"},
  {"VSCALE",		"double",	"0.001",	"Vertical Scale length for basis expansion"},
  {"TINIT",		"double",	"0.0",		"Initial time for basis"},
  {"TIME",		"double",	"0.0",		"Desired time slice"},
  {"LMAX",		"int",		"36",		"Maximum harmonic order for spherical expansion"},
  {"NMAX",		"int",		"12",		"Maximum radial order for spherical expansion"},
  {"MMAX",		"int",		"2",		"Maximum harmonic order"},
  {"NORDER",		"int",		"8",		"Number of basis functions per azimuthal harmonic"},
  {"OUTR",		"int",		"40",		"Number of radial points for output"},
  {"OUTZ",		"int",		"40",		"Number of vertical points for output"},
  {"SURFACE",		"bool",		"true",		"Make equitorial and vertical slices"},
  {"VOLUME",		"bool",		"false",	"Make volume for VTK"},
  {"AXIHGHT",		"bool",		"false",	"Compute midplane height profiles"},
  {"VHEIGHT",		"bool",		"false",	"Compute height profiles"},
  {"ALL",		"bool",		"false",	"Compute output for every time slice"},
  {"INITFLAG",		"int",		"1",		"Train set on Component (1=stars)"},
  {"PARTFLAG",		"int",		"1",		"Wakes using Component(s) [1=stars | 2=gas]"},
  {"OUTFILE",		"string",	"diskprof",	"Filename prefix"},
  {"INITIAL",		"string",	"OUT.0",	"Initial phase space file"},
  {"INFILE",		"string",	"OUT",		"Phase space file"},
  {"INDEX",		"string",	"frame.indx",	"File containing desirecd indices for PSP output"},
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
  
enum ComponentType {Star=1, Gas=2};

void add_particles(ifstream* in, PSPDump* psp, int& nbods, vector<Particle>& p)
{
  if (myid==0) {

    int nbody = nbods/numprocs;
    int nbody0 = nbods - nbody*(numprocs-1);

    MPI_Bcast(&nbody, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<Particle> t(nbody);
    vector<double> val(nbody);
    SParticle *part = psp->GetParticle(in);
    Particle bod;
    
    p = vector<Particle>(nbody0);

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
    // Send the rest of the particles
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
      MPI_Send(&val[0], nbody, MPI_DOUBLE, 1, 11, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, 1, 12, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, 1, 13, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].pos[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, 1, 14, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[0];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, 1, 15, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[1];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, 1, 16, MPI_COMM_WORLD);

      for (int i=0; i<nbody; i++) val[i] = t[i].vel[2];
      MPI_Send(&val[0], nbody, MPI_DOUBLE, 1, 17, MPI_COMM_WORLD);
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

}


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

  return -config.get<double>("ZMAX")+dz*(ipeak + delta);

}



Vector get_quart(Vector& vv, double dz)
{
  double ZMAX = config.get<double>("ZMAX");
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
  double ZMAX = config.get<double>("ZMAX");
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


void write_output(EmpCylSL& ortho, int icnt, double time)
{
  ofstream *out0, **out1;
  ofstream **out2, **out3;
  unsigned ncnt = 0;
  int nout;
  
  string OUTFILE  = config.get<string>("OUTFILE");
  double RMAX     = config.get<double>("RMAX");
  double ZMAX     = config.get<double>("ZMAX");
  int OUTR        = config.get<int   >("OUTR");
  int OUTZ        = config.get<int   >("OUTZ");
  
  bool AXIHGT     = config.get<bool>("AXIHGT");
  bool ALL        = config.get<bool>("ALL");
  bool VHEIGHT    = config.get<bool>("VHEIGHT");
  bool VOLUME     = config.get<bool>("VOLUME");
  bool SURFACE    = config.get<bool>("SURFACE");
  
  if (myid==0) {
    
    ostringstream sstr;
    sstr << "." << icnt;
    
    // ==================================================
    // Setup output files
    // ==================================================
    
    string OUTF;
    
    if (AXIHGT) {
      
      OUTF = OUTFILE + ".profile";
      if (ALL) OUTF += sstr.str();
      out0 = new ofstream(OUTF.c_str(), ios::out | ios::app);
      if (!out0) {
	cerr << "Error opening <" << OUTF << "> for output\n";
	exit(-1);
      }
      out0->setf(ios::scientific);
      out0->precision(5);
    }
    
    if (VHEIGHT) {
      
      string nsuffx[3] = {".surf", ".height", ".mid"};
      out1 = new ofstream* [3];
      
      for (int i=0; i<3; i++) {
	OUTF = OUTFILE + nsuffx[i];
	if (ALL) OUTF += sstr.str();
	out1[i] = new ofstream(OUTF.c_str());
	if (!out1[i]) {
	  cerr << "Error opening <" << OUTF << "> for output\n";
	  exit(-1);
	}
      }
      
    }
    
    nout = 5;
    string suffix[5] = {"p", "fr", "fz", "fp", "d"};
    
    if (VOLUME) {
      
      out2 = new ofstream* [nout];
      
      for (int i=0; i<nout; i++) {
	OUTF = OUTFILE + "." + suffix[i] + ".vol";
	if (ALL) OUTF += sstr.str();
	out2[i] = new ofstream(OUTF.c_str());
	if (!out2[i]) {
	  cerr << "Error opening <" << OUTF << "> for output\n";
	  exit(-1);
	}
      }
    }
    
    if (SURFACE) {
      
      out3 = new ofstream* [nout];
      
      for (int i=0; i<nout; i++) {
	OUTF = OUTFILE + "." + suffix[i] + ".surf";
	if (ALL) OUTF += sstr.str();
	out3[i] = new ofstream(OUTF.c_str());
	if (!out3[i]) {
	  cerr << "Error opening <" << OUTF << "> for output\n";
	  exit(-1);
	}
      }
    }
  }
  
  
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
      
      if ( (ncnt++)%numprocs == 0 ) {
	
	r = dR*l;
	
	for (int k=0; k<OUTZ; k++) {
	  z = -ZMAX + dz*k;
	  phi = 0.0;
	  vv[k] = ortho.accumulated_dens_eval(r, z, phi, d0);
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
      for (int l=0; l<OUTR; l++) {
	r = dR*l;
	*out0
	  << setw(15) << time
	  << setw(15) << r
	  << setw(15) << otdat[0*OUTR+l-1]
	  << setw(15) << otdat[1*OUTR+l-1]
	  << setw(15) << otdat[2*OUTR+l-1]
	  << endl;
      }
      *out0 << endl;
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
    vector<double> indat(3*OUTR*OUTR, 0.0), otdat(3*OUTR*OUTR);
    
    for (int l=0; l<OUTR; l++) {
      
      y = -RMAX + dR*l;
      
      for (int j=0; j<OUTR; j++) {
	
	x = -RMAX + dR*j;
	
	if ( (ncnt++)%numprocs == 0 ) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  for (int k=0; k<OUTZ; k++)
	    vv[k+1] = ortho.accumulated_dens_eval(r, -ZMAX + dz*k, phi, d0);
	  
	  // q = get_quart(vv, dz);
	  q = get_quart_truncated(vv, dz);
	  
	  indat[0*OUTR*OUTR + l*OUTR + j] =  get_max_dens(vv, dz);
	  indat[1*OUTR*OUTR + l*OUTR + j] =  q[1] - q[-1];
	  indat[2*OUTR*OUTR + l*OUTR + j] =  q[0];
	}
      }
      
      
      MPI_Reduce(&indat[0], &otdat[0], 3*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 0, 
		 MPI_COMM_WORLD);
      
      if (myid==0) {
	for (int i=0; i<3; i++) {
	  out1[i]->write((char *)&OUTR, sizeof(int));
	  out1[i]->write((char *)&OUTR, sizeof(int));
	}
	
	for (int l=0; l<OUTR; l++) {
	  y = -RMAX + dR*l;
	  
	  for (int j=0; j<OUTR; j++) {
	    
	    x = -RMAX + dR*j;
	    
	    if ( (ncnt++)%numprocs == 0 ) {
	      
	      r = sqrt(x*x + y*y);
	      phi = atan2(y, x);
	      
	      for (int i=0; i<3; i++) {
		out1[i]->write((char *)&x, sizeof(double));
		out1[i]->write((char *)&y, sizeof(double));
		z = otdat[i*OUTR*OUTR + l*OUTR + j];
		out1[i]->write((char *)&z, sizeof(double));
	      }
	    }
	  }
	}
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
      double p0, p, fr, fz, fp;
      
      vector<double> indat(nout*OUTZ*OUTR*OUTR, 0.0), otdat(nout*OUTZ*OUTR*OUTR);
      
      for (int k=0; k<OUTZ; k++) {
	
	z = -ZMAX + dz*k;
	
	for (int l=0; l<OUTR; l++) {
	  
	  y = -RMAX + dR*l;
	  
	  for (int j=0; j<OUTR; j++) {
	    
	    x = -RMAX + dR*j;
	    
	    r = sqrt(x*x + y*y);
	    phi = atan2(y, x);
	    
	    ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	    v = ortho.accumulated_dens_eval(r, z, phi, d0);
	    
	    indat[0*OUTZ*OUTR*OUTR + k*OUTR*OUTR + l*OUTR + j] = p;
	    indat[1*OUTZ*OUTR*OUTR + k*OUTR*OUTR + l*OUTR + j] = fr;
	    indat[2*OUTZ*OUTR*OUTR + k*OUTR*OUTR + l*OUTR + j] = fz;
	    indat[3*OUTZ*OUTR*OUTR + k*OUTR*OUTR + l*OUTR + j] = fp;
	    indat[4*OUTZ*OUTR*OUTR + k*OUTR*OUTR + l*OUTR + j] = v;
	  }
	}
      }
      
      
      MPI_Reduce(&indat[0], &otdat[0], nout*OUTZ*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
		 0, MPI_COMM_WORLD);
      
      for (int i=0; i<nout; i++) {
	out2[i]->write((char *)&OUTR, sizeof(int));
	out2[i]->write((char *)&OUTR, sizeof(int));
	out2[i]->write((char *)&OUTZ, sizeof(int));
      }
      
      int t, tn=OUTZ*OUTR*OUTR;
      
      for (int k=0; k<OUTZ; k++) {
	
	z = -ZMAX + dz*k;
	
	for (int l=0; l<OUTR; l++) {
	  
	  y = -RMAX + dR*l;
	  
	  for (int j=0; j<OUTR; j++) {
	    
	    x = -RMAX + dR*j;
	    
	    for (int i=0; i<nout; i++) {
	      out2[i]->write((char *)&x, sizeof(double));
	      out2[i]->write((char *)&y, sizeof(double));
	      out2[i]->write((char *)&z, sizeof(double));
	    }
	    
	    t = k*OUTR*OUTR + l*OUTR + j;
	    
	    out2[0]->write((char *)&(otdat[0*tn+t]) , sizeof(double));
	    out2[0]->write((char *)&valid, sizeof(int));
	    out2[1]->write((char *)&(otdat[1*tn+t]), sizeof(double));
	    out2[1]->write((char *)&valid, sizeof(int));
	    out2[2]->write((char *)&(otdat[2*tn+t]), sizeof(double));
	    out2[2]->write((char *)&valid, sizeof(int));
	    out2[3]->write((char *)&(otdat[3*tn+t]), sizeof(double));
	    out2[3]->write((char *)&valid, sizeof(int));
	    out2[4]->write((char *)&(otdat[4*tn+t]), sizeof(double));
	    out2[4]->write((char *)&valid, sizeof(int));
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
      double x, y, z=0.0, r, phi;
      double p0, d0, p, fr, fz, fp;
      
      vector<double> indat(nout*OUTR*OUTR, 0.0), otdat(nout*OUTR*OUTR);
      
      for (int l=0; l<OUTR; l++) {
	
	y = -RMAX + dR*l;
	
	for (int j=0; j<OUTR; j++) {
	  
	  if ((ncnt++)%numprocs == 0) {
	    
	    x = -RMAX + dR*j;
	    
	    r = sqrt(x*x + y*y);
	    phi = atan2(y, x);
	    
	    ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	    v = ortho.accumulated_dens_eval(r, z, phi, d0);
	    
	    indat[0*OUTR*OUTR+l*OUTR+j] = p;
	    indat[1*OUTR*OUTR+l*OUTR+j] = fr;
	    indat[2*OUTR*OUTR+l*OUTR+j] = fz;
	    indat[3*OUTR*OUTR+l*OUTR+j] = fp;
	    indat[4*OUTR*OUTR+l*OUTR+j] = v;
	  }
	}
      }
      
      MPI_Reduce(&indat[0], &otdat[0], nout*OUTZ*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
		 0, MPI_COMM_WORLD);
      
      
      if (myid==0) {
	
	for (int i=0; i<nout; i++) {
	  out3[i]->write((char *)&OUTR, sizeof(int));
	  out3[i]->write((char *)&OUTR, sizeof(int));
	  out3[i]->write((char *)&(f=-RMAX), sizeof(float));
	  out3[i]->write((char *)&(f= RMAX), sizeof(float));
	  out3[i]->write((char *)&(f=-RMAX), sizeof(float));
	  out3[i]->write((char *)&(f= RMAX), sizeof(float));
	}
	
	
	int t, tn=OUTR*OUTR;
	
	for (int l=0; l<OUTR; l++) {
	  
	  y = -RMAX + dR*l;
	  
	  for (int j=0; j<OUTR; j++) {
	    
	    x = -RMAX + dR*j;
	    
	    t = l*OUTR + j;
	    
	    out3[0]->write((char *)&(f=otdat[0*tn+t]), sizeof(float));
	    out3[1]->write((char *)&(f=otdat[1*tn+t]), sizeof(float));
	    out3[2]->write((char *)&(f=otdat[2*tn+t]), sizeof(float));
	    out3[3]->write((char *)&(f=otdat[3*tn+t]), sizeof(float));
	    out3[4]->write((char *)&(f=otdat[4*tn+t]), sizeof(float));
	  }
	}
      }
    }
    
    
    if (AXIHGT) {
      out0->close();
      delete out0;
    }
    
    if (VHEIGHT) {
      for (int i=0; i<3; i++) {
	out1[i]->close();
	delete out1[i];
      }
      delete [] out1;
    }
    
    if (VOLUME) {
      for (int i=0; i<nout; i++) {
	out2[i]->close();
	delete out2[i];
      }
      delete [] out2;
    }
    
    if (SURFACE) {
      for (int i=0; i<nout; i++) {
	out3[i]->close();
	delete out3[i];
      }
      delete [] out3;
    }

  }

}

				// Global variables
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
  // PSP input stream
  // ==================================================

  int iok = 1;
  ifstream in0, in1;
  if (myid==0) {
    in0.open(config.get<string>("INITIAL").c_str());
    if (!in0) {
      cerr << "Error opening <" << config.get<string>("INITIAL")
	   << ">" << endl;
      iok = 0;
    }

    in1.open(config.get<string>("INFILE").c_str());
    if (!in1) {
      cerr << "Error opening <" << config.get<string>("INFILE")
	   << ">" << endl;
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

  EmpCylSL::RMIN = config.get<double>("RCYLMIN");
  EmpCylSL::RMAX = config.get<double>("RCYLMAX");
  EmpCylSL::NUMX = config.get<int>("NUMX");
  EmpCylSL::NUMY = config.get<int>("NUMY");
  EmpCylSL::CMAP = true;
  EmpCylSL::logarithmic = true;
  EmpCylSL::DENS = config.get<bool>("DENS");

                                // Create expansion
				//
  EmpCylSL ortho(config.get<int>("NMAX"), 
		 config.get<int>("LMAX"), 
		 config.get<int>("MMAX"), 
		 config.get<int>("NORDER"),
		 config.get<double>("EXPSCALE"),
		 config.get<double>("VSCALE"));

  vector<Particle> particles;
  PSPDump *psp = 0;

  if (!ortho.read_cache()) {

    if (myid==0) {
      cout << "Beginning disk partition . . . " << flush;
      psp = new PSPDump (&in0, true);
      cout << "Time=" << psp->SetTime(config.get<double>("TINIT"));
    }

    if (in0.rdstate() & ios::eofbit) {
      in0.close();
      in0.open(config.get<string>("INITIAL").c_str());
    }

    partition(&in0, psp, config.get<int>("INITFLAG"), particles);
    if (myid==0) cout << "done" << endl;

    if (myid==0) cout << "Beginning disk accumulation . . . " << flush;
    ortho.setup_eof();
    ortho.setup_accumulation();
    ortho.accumulate_eof(particles);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
  
    if (myid==0) cout << "Making the EOF . . . " << flush;
    ortho.make_eof();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    if (myid==0) cout << "Making disk coefficients . . . " << flush;
    ortho.make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
  }
      

  // ==================================================
  // Open frame list
  // ==================================================

  ofstream indx(config.get<string>("INDEX").c_str());
  if (!indx) {
    cerr << "Error opening <" << config.get<string>("INDEX") 
	 << "> for output . . . continuing without writing index file" 
	 << endl;
  }

  delete psp;
  psp = new PSPDump (&in1, true);

  Dump *dump = psp->GetDump();
  bool ALL = config.get<bool>("ALL");

  if (!ALL) {
    cout << "Time=" << psp->SetTime(config.get<double>("TIME"));
    dump = psp->CurrentDump();
  }

  int icnt = 0;

  while(dump) {

    if (myid==0) {
      cout << "Time=" << dump->header.time << endl;
      cout << "Beginning disk partition . . . " << flush;
      if (ALL)
	indx << setw(15) << icnt << setw(15) << dump->header.time << endl;
    }

    if (in1.rdstate() & ios::eofbit) {
      in1.close();
      in1.open(config.get<string>("INFILE").c_str());
    }

    partition(&in1, psp, config.get<int>("PARTFLAG"), particles);

    ortho.setup_accumulation();
    ortho.accumulate(particles);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
  
    if (myid==0) cout << "Making disk coefficients . . . " << flush;
    ortho.make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;

    if (myid==0) cout << "Writing output . . . " << flush;
    write_output(ortho, icnt++, dump->header.time);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;

    dump = psp->NextDump();
  }

  MPI_Finalize();

  return 0;
}

