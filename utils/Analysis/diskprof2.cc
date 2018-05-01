/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute VTK slices, 
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
 *  MDW 11/28/08, 02/09/18
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

const std::string overview = "Compute disk potential, force and density profiles from\nPSP phase-space output files\n";

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
  
static  bool AXIHGT;
static  bool VHEIGHT;
static  bool VOLUME;
static  bool SURFACE;
static  bool VSLICE;

class Histogram
{
public:

  std::vector<double> data;
  int N;
  double R, dR, rmax;
  
  Histogram(int N, double R) : N(N), R(R)
  {
    dR = 2.0*R/(N+1);
    data.resize(N*N, 0.0);
    rmax = R + 0.5*dR;
  }

  void Reset() { std::fill(data.begin(), data.end(), 0.0); }

  void Syncr() { 
    if (myid==0)
      MPI_Reduce(MPI_IN_PLACE, &data[0], data.size(), MPI_DOUBLE, MPI_SUM, 0,
		 MPI_COMM_WORLD);
    else
      MPI_Reduce(&data[0],     &data[0], data.size(), MPI_DOUBLE, MPI_SUM, 0,
		 MPI_COMM_WORLD);
  }

  void Add(double x, double y, double m)
  {
    if (x < -rmax or x >= rmax or
	y < -rmax or y >= rmax  ) return;

    int indX = static_cast<int>(floor((x + rmax)/dR));
    int indY = static_cast<int>(floor((y + rmax)/dR));

    indX = std::min<int>(indX, N-1);
    indY = std::min<int>(indY, N-1);

    data[indY*N + indX] += m;
  }
};

enum ComponentType {Star=1, Gas=2};

void add_particles(ifstream* in, PSPDumpPtr psp, int& nbods, vector<Particle>& p,
		   Histogram& h)
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

      // Add to histogram
      //
      h.Add(part->pos(0), part->pos(1), part->mass());
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

  // Synchronize histogram
  //
  h.Syncr();
}

void partition(ifstream* in, PSPDumpPtr psp, int cflag, vector<Particle>& p,
	       Histogram& h)
{
  p.erase(p.begin(), p.end());

  int nbods = 0;
  if (cflag & Star) {
    if (myid==0) {
      nbods = psp->CurrentDump()->nstar;
      psp->GetStar();

      add_particles(in, psp, nbods, p, h);
    } else {
      add_particles(in, psp, nbods, p, h);
    }      
  }

  if (cflag & Gas) {
    if (myid==0) {
      int nbods = psp->CurrentDump()->ngas;
      psp->GetGas();

      add_particles(in, psp, nbods, p, h);
    } else {
      add_particles(in, psp, nbods, p, h);
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

  return -ZMAX + dz*(ipeak + delta);

}



Vector get_quart(Vector& vv, double dz)
{
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


void write_output(EmpCylSL& ortho, int icnt, double time, Histogram& histo)
{
  unsigned ncnt = 0;
  int nout;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << std::setfill('0') << std::setw(5) << icnt;

  nout = 7;
  string suffix[7] = {"p0", "p", "fr", "fz", "fp", "d0", "d"};

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
      
      if ( (ncnt++)%numprocs == myid ) {
	
	r = dR*l;
	
	for (int k=0; k<OUTZ; k++) {
	  z = -ZMAX + dz*k;
	  phi = 0.0;
	  vv[k+1] = ortho.accumulated_dens_eval(r, z, phi, d0);
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

      std::string OUTF = runtag + "_" + outid + "_profile" + sstr.str();
      std::ofstream out(OUTF.c_str(), ios::out | ios::app);
      if (!out) {
	cerr << "Error opening <" << OUTF << "> for output\n";
	exit(-1);
      }
      out.setf(ios::scientific);
      out.precision(5);

      for (int l=0; l<OUTR; l++) {
	r = dR*l;
	out
	  << setw(15) << time
	  << setw(15) << r
	  << setw(15) << otdat[0*OUTR+l]
	  << setw(15) << otdat[1*OUTR+l]
	  << setw(15) << otdat[2*OUTR+l]
	  << endl;
      }
      out << endl;
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
    std::vector<double> indat(3*OUTR*OUTR, 0.0), otdat(3*OUTR*OUTR);
    
    for (int j=0; j<OUTR; j++) {
	
      x = -RMAX + dR*j;
	
      for (int l=0; l<OUTR; l++) {
      
	y = -RMAX + dR*l;
      
	if ( (ncnt++)%numprocs == myid ) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  for (int k=0; k<OUTZ; k++)
	    vv[k+1] = ortho.accumulated_dens_eval(r, -ZMAX + dz*k, phi, d0);
	  
	  // q = get_quart(vv, dz);
	  q = get_quart_truncated(vv, dz);
	  
	  indat[(0*OUTR + j)*OUTR + l] =  get_max_dens(vv, dz);
	  indat[(1*OUTR + j)*OUTR + l] =  q[1] - q[-1];
	  indat[(2*OUTR + j)*OUTR + l] =  q[0];
	}
      }
    }
      
      
    MPI_Reduce(&indat[0], &otdat[0], 3*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
      
    if (myid==0) {

      VtkGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::string names[3] = {"surf", "height", "mid"};
      
      std::vector<double> data(OUTR*OUTR);

      for (int i=0; i<3; i++) {

	for (int l=0; l<OUTR; l++) {
	  for (int j=0; j<OUTR; j++) {
	    data[i*OUTR + j] = otdat[(i*OUTR + l)*OUTR + j];
	  }
	}

	vtk.Add(data, names[i]);
      }

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_posn" + sstr.str();
      vtk.Write(sout.str());
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
    double p0, d0, p, fr, fz, fp;
    
    size_t blSiz = OUTZ*OUTR*OUTR;
    vector<double> indat(nout*blSiz, 0.0), otdat(nout*blSiz);
    
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
    
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTZ*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      VtkGrid vtk(OUTR, OUTR, OUTZ, -RMAX, RMAX, -RMAX, RMAX, -ZMAX, ZMAX);

      std::vector<double> data(OUTR*OUTR*OUTZ);

      for (int n=0; n<nout; n++) {

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
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      VtkGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::vector<double> data(OUTR*OUTR);

      for (int n=0; n<nout; n++) {
	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    data[l*OUTR + j] = otdat[(n*OUTR+j)*OUTR+l];
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      vtk.Add(histo.data, "histo");

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_surface" + sstr.str();
      vtk.Write(sout.str());
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
    
    std::vector<double> indat(nout*OUTR*OUTR, 0.0), otdat(nout*OUTR*OUTR);
    
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
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTZ,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      VtkGrid vtk(OUTR, OUTZ, 1, -RMAX, RMAX, -ZMAX, ZMAX, 0, 0);

      std::vector<double> data(OUTR*OUTZ);

      for (int n=0; n<nout; n++) {
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


int
main(int argc, char **argv)
{
  int nice, numx, numy, lmax, mmax, nmax, norder;
  int initc, partc, beg, end, stride, init;
  double rcylmin, rcylmax, rscale, vscale;
  bool DENS, HEIGHT, PCA, PVD, verbose = false, mask = false;
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
    ("surface",
     po::value<bool>(&SURFACE)->default_value(true),
     "make equitorial slices")
    ("vslice",
     po::value<bool>(&VSLICE)->default_value(true),
     "make vertical slices")
    ("volume",
     po::value<bool>(&VOLUME)->default_value(false),
     "make volume for VTK rendering")
    ("axihgt",
     po::value<bool>(&AXIHGT)->default_value(false),
     "compute midplane height profiles")
    ("height",
     po::value<bool>(&HEIGHT)->default_value(false),
     "compute height profiles")
    ("pca",
     po::value<bool>(&PCA)->default_value(false),
     "perform the PCA analysis for the disk")
    ("density",
     po::value<bool>(&DENS)->default_value(true),
     "compute density")
    ("pvd",
     po::value<bool>(&PVD)->default_value(false),
     "Compute PVD file for ParaView")
    ("initcomp",
     po::value<int>(&initc)->default_value(1),
     "train on Component (1=stars)")
    ("partflag",
     po::value<int>(&partc)->default_value(1),
     "Wakes using Component(s) [1=stars | 2=gas]")
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
    
  if (PCA) {
    EmpCylSL::SELECT = true;
      ortho.setHall(runtag + "_" + outid + ".pca", 1);
  }

  vector<Particle> particles;
  PSPDumpPtr psp;
  Histogram histo(OUTR, RMAX);
  
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
      
    partition(&in0, psp, initc, particles, histo);
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
    
    int icnt = 0;
    
    //------------------------------------------------------------ 

    if (myid==0) {
      tnow = dump->header.time;
      cout << "Beginning disk partition [time=" << tnow
	   << ", index=" << indx << "] . . . "  << flush;
      times.push_back(tnow);
      ostringstream filen;
      filen << runtag << "_" << outid << "_surface."
	    << std::setfill('0') << std::setw(5) << indx << ".vtr";
      outfiles.push_back(filen.str());
    }
      
    if (in1.rdstate() & ios::eofbit) {
      in1.close();
      in1.open(s1.str());
    }
      
    histo.Reset();		// Reset surface histogram

    partition(&in1, psp, partc, particles, histo);
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
      
    if (myid==0) cout << "Writing output . . . " << flush;
    write_output(ortho, indx, dump->header.time, histo);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
    
  } // Dump loop

  // Create PVD file
  //
  if (myid==0 and PVD) {
    writePVD(runtag+".pvd", times, outfiles);
  }

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

