/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute VTK slices, and compute volume
 *  for VTK rendering
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
 *  MDW 11/28/08, 02/09/18, 11/21/19
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
#include <PSP2.H>
#include <interp.h>
#include <EmpCylSL.h>

#include <localmpi.h>
#include <foarray.H>

#include <VtkGrid.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

const std::string overview = "Compute disk potential, force and density profiles from\nPSP phase-space output files";

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
static  bool PROBE;

std::vector<double> c0 = {0.0, 0.0, 0.0};

// Temporary center offset
// std::vector<double> c0 = {0.0, 0.0, -0.00217};

class Histogram
{
public:

  std::vector<double> dataXY, dataXZ, dataYZ;
  std::vector<std::vector<double>> dataZ;
  int N, M;
  double R, dR, rmax;
  double Z, dZ, zmax;
  
  Histogram(int N, int M, double R, double Z) : N(N), M(M)
  {
    dR = 2.0*R/(N+1);
    dZ = 2.0*Z/(M+1);

    rmax = R + 0.5*dR;
    zmax = Z + 0.5*dZ;

    dataXY.resize(N*N, 0.0);
    dataXZ.resize(N*M, 0.0);
    dataYZ.resize(N*M, 0.0);
    dataZ .resize(N*N);

  }

  void Reset() {
    std::fill(dataXY.begin(), dataXY.end(), 0.0);
    std::fill(dataXZ.begin(), dataXZ.end(), 0.0);
    std::fill(dataYZ.begin(), dataYZ.end(), 0.0);
    for (auto & v : dataZ) v.clear();
  }

  void Syncr() { 
    if (myid==0) {
      MPI_Reduce(MPI_IN_PLACE, &dataXY[0], dataXY.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &dataXZ[0], dataXZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &dataYZ[0], dataYZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else {
      MPI_Reduce(&dataXY[0],   &dataXY[0], dataXY.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&dataXZ[0],   &dataXZ[0], dataXZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&dataYZ[0],   &dataYZ[0], dataXZ.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    std::vector<double> work;
    int nz;

    for (int n=1; n<numprocs; n++) {
      if (myid==0) {
	for (auto & v : dataZ) {
	  MPI_Recv(&nz, 1, MPI_INT, n, 110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  work.resize(nz);
	  MPI_Recv(&work[0], nz, MPI_DOUBLE, n, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  v.insert(std::end(v), std::begin(work), std::end(work));
	}
      } else if (myid==n) {
	for (auto v : dataZ) {
	  nz = v.size();
	  MPI_Send(&nz, 1,MPI_INT, 0, 110, MPI_COMM_WORLD);
	  MPI_Send(&v[0], nz, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
	}
      }
    }

    if (myid==0) {
      for (auto & v : dataZ) std::sort(std::begin(v), std::end(v));
    }
  }

  void Add(double x, double y, double z, double m)
  {
    if (x < -rmax or x >= rmax or
	y < -rmax or y >= rmax or
	z < -zmax or z >= zmax) return;

    int indX = static_cast<int>(floor((x + rmax)/dR));
    int indY = static_cast<int>(floor((y + rmax)/dR));
    int indZ = static_cast<int>(floor((z + zmax)/dZ));

    indX = std::min<int>(indX, N-1);
    indY = std::min<int>(indY, N-1);
    indZ = std::min<int>(indZ, M-1);

    dataXY[indY*N + indX] += m;
    dataXZ[indX*M + indZ] += m;
    dataYZ[indY*M + indZ] += m;
    dataZ [indY*N + indX].push_back(z);
  }
};

void add_particles(PSPptr psp, int& nbods, vector<Particle>& p, Histogram& h)
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

      // Make a Particle
      //
      bod.mass = part->mass();
      for (int k=0; k<3; k++) bod.pos[k] = part->pos(k) - c0[k];
      for (int k=0; k<3; k++) bod.vel[k] = part->vel(k);
      bod.indx = part->indx();
      p.push_back(bod);

      part = psp->NextParticle();

      // Add to histogram
      //
      if (part) h.Add(part->pos(0) - c0[0],
		      part->pos(1) - c0[1],
		      part->pos(2) - c0[2],
		      part->mass());
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
	for (int k=0; k<3; k++) t[i].pos[k] = part->pos(k) - c0[k];
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

    // Add to histogram
    //
    for (auto part : t)
      h.Add(part.pos[0], part.pos[1], part.pos[2], part.mass);
  }

  // Synchronize histogram
  //
  h.Syncr();
}

void partition(PSPptr psp, std::string& name, vector<Particle>& p, Histogram& h)
{
  p.erase(p.begin(), p.end());

  PSPstanza* stanza;

  int iok = 1, nbods = 0;

  if (myid==0) {
    stanza = psp->GetNamed(name);
    if (stanza==0)
      iok = 0;
    else
      nbods = stanza->comp.nbod;
  }

  MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (iok==0) {
    if (myid==0) std::cerr << "Could not find component named <" << name << ">" << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  add_particles(psp, nbods, p, h);
}

// Find peak density
// -----------------
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
  double del  = vv[ipeak+1] - vv[ipeak-1];
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
  int noutV = 7, noutS = 10;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << std::setfill('0') << std::setw(5) << icnt;

  string suffix[10] = {"p0", "p", "fr", "fz", "fp", "d0", "d",
		       "z10", "z50", "z90"};

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
	    data[l*OUTR + j] = otdat[(i*OUTR + l)*OUTR + j];
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

      VtkGrid vtk(OUTR, OUTR, OUTZ, -RMAX, RMAX, -RMAX, RMAX, -ZMAX, ZMAX);

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
      
      for (int j=0; j<OUTR; j++) {
	for (int l=0; l<OUTR; l++) {
	  otdat[(7*OUTR+j)*OUTR+l] = 0.0;
	  otdat[(8*OUTR+j)*OUTR+l] = 0.0;
	  otdat[(9*OUTR+j)*OUTR+l] = 0.0;
	  
	  // Check for number in the histogram bin
	  //
	  int numZ = histo.dataZ[l*OUTR+j].size();
	  if (numZ>0) {
	    otdat[(7*OUTR+j)*OUTR+l] = histo.dataZ[l*OUTR+j][floor(0.1*numZ)];
	    otdat[(8*OUTR+j)*OUTR+l] = histo.dataZ[l*OUTR+j][floor(0.5*numZ)];
	    otdat[(9*OUTR+j)*OUTR+l] = histo.dataZ[l*OUTR+j][floor(0.9*numZ)];
	  }
	}
      }

      VtkGrid vtkXY(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::vector<double> dataXY(OUTR*OUTR);

      for (int n=0; n<noutS; n++) {
	for (int j=0; j<OUTR; j++) {
	  for (int l=0; l<OUTR; l++) {
	    dataXY[j*OUTR + l] = otdat[(n*OUTR+j)*OUTR+l];
	  }
	}
	vtkXY.Add(dataXY, suffix[n]);
      }

      vtkXY.Add(histo.dataXY, "histoXY");

      std::ostringstream sout;
      sout << runtag + "_" + outid + "_surface" + sstr.str();
      vtkXY.Write(sout.str());
      
      VtkGrid vtkZ(OUTR, OUTZ, 1, -RMAX, RMAX, -ZMAX, ZMAX, 0, 0);

      vtkZ.Add(histo.dataXZ, "histoXZ");
      vtkZ.Add(histo.dataYZ, "histoYZ");

      sout.str("");
      sout << runtag + "_" + outid + "_vsurf" + sstr.str();
      vtkZ.Write(sout.str());
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
      
      VtkGrid vtk(OUTR, OUTZ, 1, -RMAX, RMAX, -ZMAX, ZMAX, 0, 0);

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
    }
  } // END: VSLICE


  if (PROBE) {
    
    // ==================================================
    // Write 1d prob parallel and perpendicular to
    // to disk plane
    // ==================================================
    
    double v;
    float f;
    
    double dR = RMAX/OUTR;
    double z=0, r=0, phi=0;
    double p0, d0, p, fr, fz, fp;
    
    std::vector<double> indat(3*OUTR, 0.0), otdat(3*OUTR);
    
    for (int j=1; j<=OUTR; j++) {
      r = dR*j;

      if ((ncnt++)%numprocs == myid) {
	  
	phi = 0.0;
	ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	
	indat[0*OUTR + j] = fr;
	
	phi = 0.5*M_PI;
	ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	indat[1*OUTR + j] = fr;
	
	phi = 0.0;
	ortho.accumulated_eval(z, r, phi, p0, p, fr, fz, fp);
	
	indat[2*OUTR + j] = fz;
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], 3*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      std::ostringstream sout;
      sout << runtag + "_" + outid + ".probe";
      std::ofstream out(sout.str());
      if (out) {
	// Begin: header
	//
	out << std::right
	    << "# " << std::setw(14) << "r |"
	    << std::setw(16) << "f_r(phi=0) |"
	    << std::setw(16) << "f_r(phi=pi/2) |"
	    << std::setw(16) << "f_z(r=0) |"
	    << std::endl;

	out << std::right
	    << "# " << std::setw(14) << "[1] |"
	    << std::setw(16) << "[2] |"
	    << std::setw(16) << "[3] |"
	    << std::setw(16) << "[4] |"
	    << std::endl;

	out << std::right
	    << "#" << std::setfill('-') << std::setw(15) << "+"
	    << std::setw(16) << "+"
	    << std::setw(16) << "+"
	    << std::setw(16) << "+"
	    << std::endl << std::setfill(' ');
	//
	// END: header

	for (int j=1; j<=OUTR; j++) {
	  r = dR*j;
	  out << std::setw(16) << r
	      << std::setw(16) << otdat[0*OUTR + j]
	      << std::setw(16) << otdat[1*OUTR + j]
	      << std::setw(16) << otdat[2*OUTR + j]
	      << std::endl;
	}
      }
    }
  } // END: PROBE

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
  int beg, end, stride, init;
  double rcylmin, rcylmax, rscale, vscale;
  bool DENS, PCA, PVD, verbose = false, mask = false, cmap, logl, ignore;
  std::string CACHEFILE, cname, pname, dir("./");

  //
  // Parse Command line
  //
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "produce this help message")
    ("noCommand,X",
     "do not save command line")
    ("verbose,v",
     "verbose output")
    ("mask,b",
     "blank empty cells")
    ("OUT",
     "assume original, single binary PSP files as input")
    ("SPL",
     "assume new split binary PSP files as input")
    ("nice",
     po::value<int>(&nice)->default_value(0), 
     "number of bins in x direction")
    ("x,x",
     po::value<double>(&c0[0])->default_value(0.0), 
     "x-position offset for phase space")
    ("y,y",
     po::value<double>(&c0[1])->default_value(0.0), 
     "y-position offset for phase space")
    ("z,z",
     po::value<double>(&c0[2])->default_value(0.0), 
     "x-position offset for phase space")
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
     "make equatorial slices")
    ("vslice",
     po::value<bool>(&VSLICE)->default_value(true),
     "make vertical slices")
    ("probe",
     po::value<bool>(&PROBE)->default_value(true),
     "make 1d cuts in and perpendicular to the equitorial plane")
    ("volume",
     po::value<bool>(&VOLUME)->default_value(false),
     "make volume for VTK rendering")
    ("axihgt",
     po::value<bool>(&AXIHGT)->default_value(false),
     "compute midplane height profiles")
    ("height",
     po::value<bool>(&VHEIGHT)->default_value(false),
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
    ("compname",
     po::value<std::string>(&cname)->default_value("stars"),
     "train on Component (default=stars)")
    ("partname",
     po::value<std::string>(&pname)->default_value("stars"),
     "Wakes using Component name")
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
     po::value<std::string>(&outid)->default_value("diskprof4"),
     "Filename prefix")
    ("cachefile",
     po::value<std::string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "cachefile name")
    ("runtag",
     po::value<std::string>(&runtag)->default_value("run1"),
     "runtag for phase space files")
    ("dir,d",
     po::value<std::string>(&dir),
     "directory for SPL files")
    ("cmap",
     po::value<bool>(&cmap)->default_value(true),
     "map radius into semi-infinite interval in cylindrical grid computation")
    ("ignore",
     po::value<bool>(&ignore)->default_value(false),
     "rebuild EOF grid if input parameters do not match the cachefile")
    ("logl",
     po::value<bool>(&logl)->default_value(true),
     "use logarithmic radius scale in cylindrical grid computation")
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
    std::cout << std::string(60, '-') << std::endl;
    std::cout << overview << std::endl;
    std::cout << std::string(60, '-') << std::endl << std::endl;
    std::cout << desc     << std::endl;
    return 1;
  }
 
  if (vm.count("verbose")) verbose = true;

  if (vm.count("mask")) mask = true;

  bool SPL = false;
  if (vm.count("SPL")) SPL = true;
  if (vm.count("OUT")) SPL = false;

  if (vm.count("noCommand")==0) {
    std::string cmdFile = runtag + "." + outid + ".cmd_line";
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
  std::ostringstream s0, s1;
  if (myid==0) {
    if (SPL) s0 << "SPL.";
    else     s0 << "OUT.";
    s0 << runtag << "."
       << std::setw(5) << std::setfill('0') << init;

    std::string file = dir + s0.str();
    std::ifstream in(file);
    if (!in) {
      cerr << "Error opening <" << file << ">" << endl;
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

  // Set parameters from the given CACHEFILE
  //
  if (not ignore) {

    std::ifstream in(CACHEFILE);
    if (!in) {
      std::cerr << "Error opening cachefile named <" 
		<< CACHEFILE << "> . . ."
		<< std::endl
		<< "I will build <" << CACHEFILE
		<< "> but it will take some time."
		<< std::endl
		<< "If this is NOT what you want, "
		<< "stop this routine and specify the correct file."
		<< std::endl;
    } else {

      int tmp;
    
      in.read((char *)&mmax,    sizeof(int));
      in.read((char *)&numx,    sizeof(int));
      in.read((char *)&numy,    sizeof(int));
      in.read((char *)&nmax,    sizeof(int));
      in.read((char *)&norder,  sizeof(int));
      
      in.read((char *)&tmp,     sizeof(int)); 
      if (tmp) DENS = true;
      else     DENS = false;
      
      in.read((char *)&tmp,     sizeof(int)); 
      if (tmp) cmap = true;
      else     cmap = false;

      in.read((char *)&rcylmin, sizeof(double));
      in.read((char *)&rcylmax, sizeof(double));
      in.read((char *)&rscale,  sizeof(double));
      in.read((char *)&vscale,  sizeof(double));
    }
  }

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAP        = cmap;
  EmpCylSL::logarithmic = logl;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::CACHEFILE   = CACHEFILE;

				// Create expansion
				//
  EmpCylSL ortho(nmax, lmax, mmax, norder, rscale, vscale);
    
  if (PCA) {
    EmpCylSL::PCAVAR = true;
      ortho.setHall(runtag + "_" + outid + ".pca", 1);
  }

  vector<Particle> particles;
  PSPptr psp;
  Histogram histo(OUTR, OUTZ, RMAX, ZMAX);
  
  std::vector<double> times;
  std::vector<std::string> outfiles;

  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==0) {
    
    if (myid==0) {
      if (SPL) psp = std::make_shared<PSPspl>(s0.str(), dir, true);
      else     psp = std::make_shared<PSPout>(s0.str(), true);
      std::cout << "Beginning disk partition [time="
		<< psp->CurrentTime()
		<< "] . . . " << std::flush;
    }
      
    partition(psp, cname, particles, histo);

    if (myid==0)
      std::cout << "done" << endl << endl
		<< setw(4) << "#" << "  " << setw(8) << "Number" << endl 
		<< setfill('-')
		<< setw(4) << "-" << "  " << setw(8) << "-" << endl
		<< setfill(' ');

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


  std::string file;

  for (int indx=beg; indx<=end; indx+=stride) {

    // ==================================================
    // PSP input stream
    // ==================================================

    iok = 1;
    if (myid==0) {
      s1.str("");		// Clear stringstream
      if (SPL) s1 << "SPL.";
      else     s1 << "OUT.";
      s1 << runtag << "."<< std::setw(5) << std::setfill('0') << indx;
      
				// Check for existence of next file
      file = dir + s1.str();
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
    if (myid==0) {

      if (SPL) psp = std::make_shared<PSPspl>(s1.str(), dir, true);
      else     psp = std::make_shared<PSPout>(file, true);

      tnow = psp->CurrentTime();
      cout << "Beginning disk partition [time=" << tnow
	   << ", index=" << indx << "] . . . "  << flush;
      times.push_back(tnow);
      ostringstream filen;
      filen << runtag << "_" << outid << "_surface."
	    << std::setfill('0') << std::setw(5) << indx << ".vtr";
      outfiles.push_back(filen.str());
    }
      
    histo.Reset();		// Reset surface histogram

    partition(psp, pname, particles, histo);
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
    double time = 0.0;
    if (myid==0) time = psp->CurrentTime();
    write_output(ortho, indx, time, histo);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
    
  } // Dump loop

  // Create PVD file
  //
  if (myid==0 and PVD) {
    writePVD(runtag+"." + outid + ".pvd", times, outfiles);
  }

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
}

