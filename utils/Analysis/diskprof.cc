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
#include <PSP.H>
#include <interp.H>
#include <EmpCylSL.H>

#include <localmpi.H>
#include <foarray.H>

#include <VtkGrid.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif
  
class Histogram
{
public:

  std::vector<double> data;
  int N;
  double R, dR, rmax;
  
  Histogram(int N, double R) : N(N), R(R)
  {
    N = std::max<int>(N, 2);
    dR = 2.0*R/(N-1);
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

    data[indX*N + indY] += m;
  }
};

void add_particles(PSPptr psp, vector<Particle>& p, Histogram& h)
{
  if (myid==0) {

    int nbods = psp->GetStanza()->comp.nbod;
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
      for (int k=0; k<3; k++) bod.pos[k] = part->pos(k);
      for (int k=0; k<3; k++) bod.vel[k] = part->vel(k);
      bod.indx = part->indx();
      p.push_back(bod);

      part = psp->NextParticle();

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
  }

  // Synchronize histogram
  //
  h.Syncr();
}

void partition(PSPptr psp, std::string& comp, vector<Particle>& p, Histogram& h)
{
  if (not psp->GetNamed(comp)) {
    if (myid==0) 
      std::cout << "No component named <" << comp << ">" << std::endl;
    MPI_Finalize();
    exit(-2);
  }
  p.erase(p.begin(), p.end());
  add_particles(psp, p, h);
}

double ZMAX, RMAX;
string OUTFILE;
int OUTR, OUTZ;
bool AXIHGT, ALL, VHEIGHT, VOLUME, SURFACE;

				// Find peak density

double get_max_dens(Eigen::VectorXd& vv, double dz)
{
  int sz = vv.size();

  int ipeak=0;
  double pvalue=-1.0e18;

  for (int i=0; i<sz; i++) {
    if (vv[i] > pvalue) {
      ipeak = i;
      pvalue = vv[i];
    }
  }

  if (ipeak == 0   ) ipeak = 1;
  if (ipeak == sz-1) ipeak = sz-2;

				// Solution of 2nd order Lagrange interpolation
				// formula
  double delta;
  double del =  vv[ipeak+1] - vv[ipeak-1];
  double ddel = vv[ipeak+1] + vv[ipeak-1] - 2.0*vv[ipeak];

  if (fabs(ddel) < 1.0e-4) 
    delta = 0.0;		// We're at the peak!
  else
    delta = - 0.5*del/ddel;

  return -ZMAX+dz*(ipeak + delta);

}



Eigen::VectorXd get_quart(Eigen::VectorXd& vv, double dz)
{
  int sz = vv.size();

  double next, prev=0.0;
  Eigen::VectorXd sum(sz), zz(sz);

  for (int i=0; i<sz; i++) {

    if (vv[i] > 0.0)
      next = vv[i];
    else
      next = 0.0;
    
    zz [i] = -ZMAX + dz*i;
    sum[i] = 0.5*(prev + next);	// Trapezoidal rule
    if (i>0) sum[i] += sum[i-1];
    prev = next;
  }

  double max = sum[sz-1];
  Eigen::VectorXd ret(3);
  
  ret[0] = odd2(0.25*max, sum, zz);
  ret[1] = odd2(0.50*max, sum, zz);
  ret[2] = odd2(0.75*max, sum, zz);

  return ret;
}

Eigen::VectorXd get_quart_truncated(Eigen::VectorXd& vv, double dz)
{
  int sz = vv.size();

  int ipeak=0;			// First find peak
  double pvalue=-1.0e18;

  for (int i=0; i<sz; i++) {
    if (vv[i] > pvalue) {
      ipeak = i;
      pvalue = vv[i];
    }
  }
				// Zero out above and below first zero
				// from peak

  int lo1 = ipeak;
  int hi1 = ipeak;

  for (; lo1>0; lo1--) {
    if (vv[lo1]<0.0) break;
  }

  for (; hi1<sz; hi1++) {
    if (vv[hi1]<0.0) break;
  }

  double next, prev=0.0;
  int sz1 = hi1 - lo1 + 1;
  Eigen::VectorXd sum(sz1), zz(sz1);

  for (int i=0; i<sz1; i++) {

    if (vv[i+lo1] > 0.0)
      next = vv[i+lo1];
    else
      next = 0.0;
    
    zz [i] = -ZMAX + dz*(i+lo1);
    sum[i] = 0.5*(prev + next);	// Trapezoidal rule
    if (i>lo1) sum[i] += sum[i-1];
    prev = next;
  }

  double max = sum[sz1-1];
  Eigen::VectorXd ret(3);
  
  ret[0] = odd2(0.25*max, sum, zz);
  ret[1] = odd2(0.50*max, sum, zz);
  ret[2] = odd2(0.75*max, sum, zz);

  return ret;
}


void write_output(EmpCylSL& ortho, double time, Histogram& histo)
{
  unsigned ncnt = 0;
  int nout;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;

  nout = 7;
  string suffix[7] = {"p0", "p", "fr", "fz", "fp", "d0", "d"};

  // ==================================================
  // Axisymmetric structure (for GNUPLOT)
  // ==================================================

  if (AXIHGT) {
    
    double dR = RMAX/(OUTR-1);
    double dz = ZMAX/(OUTZ-1);
    double z, r, phi, hh, d0;
    Eigen::VectorXd vv(OUTZ);
    Eigen::VectorXd q;
    
    vector<double> indat(3*OUTR, 0.0), otdat(3*OUTR);
    
    for (int l=0; l<OUTR; l++) {
      
      if ( (ncnt++)%numprocs == myid ) {
	
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

      string OUTF = OUTFILE + ".profile";
      if (ALL) OUTF += sstr.str();
      ofstream out(OUTF.c_str(), ios::out | ios::app);
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
    Eigen::VectorXd vv(OUTZ);
    Eigen::VectorXd q;
    vector<double> indat(3*OUTR*OUTR, 0.0), otdat(3*OUTR*OUTR);
    
    for (int l=0; l<OUTR; l++) {
      
      y = -RMAX + dR*l;
      
      for (int j=0; j<OUTR; j++) {
	
	x = -RMAX + dR*j;
	
	if ( (ncnt++)%numprocs == myid ) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  for (int k=0; k<OUTZ; k++)
	    vv[k] = ortho.accumulated_dens_eval(r, -ZMAX + dz*k, phi, d0);
	  
	  // q = get_quart(vv, dz);
	  q = get_quart_truncated(vv, dz);
	  
	  indat[(0*OUTR + l)*OUTR + j] =  get_max_dens(vv, dz);
	  indat[(1*OUTR + l)*OUTR + j] =  q[1] - q[-1];
	  indat[(2*OUTR + l)*OUTR + j] =  q[0];
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
      sout << OUTFILE + "_posn";
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
	  
	  indat[((0*OUTZ + k)*OUTR + l)*OUTR + j] = p0;
	  indat[((1*OUTZ + k)*OUTR + l)*OUTR + j] = p;
	  indat[((2*OUTZ + k)*OUTR + l)*OUTR + j] = fr;
	  indat[((3*OUTZ + k)*OUTR + l)*OUTR + j] = fz;
	  indat[((4*OUTZ + k)*OUTR + l)*OUTR + j] = fp;
	  indat[((5*OUTZ + k)*OUTR + l)*OUTR + j] = d0;
	  indat[((6*OUTZ + k)*OUTR + l)*OUTR + j] = v;
	}
      }
    }
    
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTZ*OUTR*OUTR, MPI_DOUBLE, MPI_SUM, 
	       0, MPI_COMM_WORLD);
    
    if (myid==0) {

      VtkGrid vtk(OUTR, OUTR, OUTZ, -RMAX, RMAX, -RMAX, RMAX, -ZMAX, ZMAX);

      std::vector<double> data(OUTR*OUTR*OUTZ);

      for (int n=0; n<nout; n++) {

	for (int k=0; k<OUTZ; k++) {
	  for (int l=0; l<OUTR; l++) {
	    for (int j=0; j<OUTR; j++) {
	      data[(j*OUTR + l)*OUTR + k] = otdat[((n*OUTZ + k)*OUTR + l)*OUTR + j];
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
    double x, y, z=0.0, r, phi;
    double p0, d0, p, fr, fz, fp;
    
    vector<double> indat(nout*OUTR*OUTR, 0.0), otdat(nout*OUTR*OUTR);
    
    for (int l=0; l<OUTR; l++) {
      
      y = -RMAX + dR*l;
      
      for (int j=0; j<OUTR; j++) {
	
	if ((ncnt++)%numprocs == myid) {
	  
	  x = -RMAX + dR*j;
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  ortho.accumulated_eval(r, z, phi, p0, p, fr, fz, fp);
	  v = ortho.accumulated_dens_eval(r, z, phi, d0);
	  
	  indat[(0*OUTR+l)*OUTR+j] = p0;
	  indat[(1*OUTR+l)*OUTR+j] = p;
	  indat[(2*OUTR+l)*OUTR+j] = fr;
	  indat[(3*OUTR+l)*OUTR+j] = fz;
	  indat[(4*OUTR+l)*OUTR+j] = fp;
	  indat[(5*OUTR+l)*OUTR+j] = d0;
	  indat[(6*OUTR+l)*OUTR+j] = v;
	}
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], nout*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      VtkGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

      std::vector<double> data(OUTR*OUTR);

      for (int n=0; n<nout; n++) {
	for (int l=0; l<OUTR; l++) {
	  for (int j=0; j<OUTR; j++) {
	    data[l*OUTR + j] = otdat[(n*OUTR+l)*OUTR+j];
	  }
	}
	vtk.Add(data, suffix[n]);
      }

      vtk.Add(histo.data, "histo");

      std::ostringstream sout;
      sout << OUTFILE + "_surface";
      vtk.Write(sout.str());
    }
  }
}


int
main(int argc, char **argv)
{
  
#ifdef DEBUG
  sleep(20);
#endif  
  
  int NICE, NUMX, NUMY, NMAX, LMAX, MMAX, NORDER, OUTR, PARTFLAG;
  double RCYLMIN, RCYLMAX, RSCALE, VSCALE, RMAX, TINIT;
  std::string INITIAL, INFILE, CACHEFILE, OUTFILE, compname;
  bool DENS, PCA;

  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  po::options_description desc("Compute disk potential, force and density profiles\nfrom PSP phase-space output files\n\nAllowed options");
  desc.add_options()
    ("help,h",                                                                       "Print this help message")
    ("OUT",
     "assume that PSP files are in original format")
    ("SPL",
     "assume that PSP files are in split format")
    ("NICE",                po::value<int>(&NICE)->default_value(0),
     "system priority")
    ("DENS",                po::value<bool>(&DENS)->default_value(true),
     "compute density")
    ("RMAX",                po::value<double>(&RMAX)->default_value(0.1),
     "maximum radius for output")
    ("ZMAX",                po::value<double>(&ZMAX)->default_value(0.1),
     "maximum height for output")
    ("RCYLMIN",             po::value<double>(&RCYLMIN)->default_value(0.001),
     "number of scale lengths for minimum radius in table")
    ("RCYLMAX",             po::value<double>(&RCYLMAX)->default_value(20.0),
     "number of scale lengths for maximum radius in table")
    ("NUMX",                po::value<int>(&NUMX)->default_value(128),
     "number of radial table entries")
    ("NUMY",                po::value<int>(&NUMY)->default_value(64),
     "number of vertical table entries")
    ("RSCALE",              po::value<double>(&RSCALE)->default_value(0.01),
     "Radial scale length for basis expansion")
    ("VSCALE",              po::value<double>(&VSCALE)->default_value(0.001),
     "Vertical Scale length for basis expansion")
    ("TINIT",               po::value<double>(&TINIT)->default_value(0.0),
     "Initial time for PSP dump")
    ("LMAX",                po::value<int>(&LMAX)->default_value(36),
     "Maximum harmonic order for spherical expansion")
    ("NMAX",                po::value<int>(&NMAX)->default_value(8),
     "Maximum radial order for spherical expansion")
    ("MMAX",                po::value<int>(&MMAX)->default_value(4),
     "Maximum harmonic order")
    ("NORDER",              po::value<int>(&NORDER)->default_value(4),
     "Number of basis functions per azimuthal harmonic")
    ("OUTR",                po::value<int>(&OUTR)->default_value(40),
     "Number of radial points for output")
    ("OUTZ",                po::value<int>(&OUTZ)->default_value(40),
     "Number of vertical points for output")
    ("SURFACE",             po::value<bool>(&SURFACE)->default_value(true),
     "Make equatorial and vertical slices")
    ("VOLUME",              po::value<bool>(&VOLUME)->default_value(false),
     "Make volume for VTK")
    ("AXIHGT",              po::value<bool>(&AXIHGT)->default_value(false),
     "Compute midplane height profiles")
    ("VHEIGHT",             po::value<bool>(&VHEIGHT)->default_value(false),
     "Compute height profiles")
    ("ALL",                 po::value<bool>(&ALL)->default_value(false),
     "Compute output for every time slice")
    ("PCA",                 po::value<bool>(&PCA)->default_value(false),
     "Perform the PCA analysis for the disk")
    ("compname",            po::value<std::string>(&compname),
     "Component name (no default)")
    ("PARTFLAG",            po::value<int>(&PARTFLAG)->default_value(1),
     "Wakes using Component(s) [1=stars | 2=gas]")
    ("OUTFILE",             po::value<string>(&OUTFILE)->default_value("diskprof"),
     "Filename prefix")
    ("CACHEFILE",           po::value<string>(&CACHEFILE)->default_value(".eof.cache.file"),
     "Cachefile name")
    ("INITIAL",             po::value<string>(&INITIAL)->default_value("OUT.0"),
     "Initial phase space file")
    ("INFILE",              po::value<string>(&INFILE)->default_value("OUT"),
     "Phase space file")
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

  if (vm.count("help")) {
    std::cout << std::string(60, '-') << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }
 

  // ==================================================
  // Nice process
  // ==================================================

  if (NICE>0)
    setpriority(PRIO_PROCESS, 0, NICE);


  // ==================================================
  // PSP input stream
  // ==================================================

  int iok = 1;
  if (myid==0) {
    std::ifstream in(INITIAL);
    if (!in) {
      cerr << "Error opening <" << INITIAL << ">" << endl;
      iok = 0;
    }
    in.close();
    in.open(INFILE);
    if (!in) {
      cerr << "Error opening <" << INFILE << ">" << endl;
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

  EmpCylSL::RMIN        = RCYLMIN;
  EmpCylSL::RMAX        = RCYLMAX;
  EmpCylSL::NUMX        = NUMX;
  EmpCylSL::NUMY        = NUMY;
  EmpCylSL::CMAPR       = 1;
  EmpCylSL::CMAPZ       = 1;
  EmpCylSL::logarithmic = true;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::CACHEFILE   = CACHEFILE;

                                // Create expansion
				//
  EmpCylSL ortho(NMAX, LMAX, MMAX, NORDER, RSCALE, VSCALE);

  if (PCA) {
    EmpCylSL::PCAVAR = true;
    ortho.setHall(OUTFILE + ".pca", 1);
  }

  vector<Particle> particles;
  Histogram histo(OUTR, RMAX);
  PSPptr psp;

  if (ortho.read_cache()==0) {

    //------------------------------------------------------------ 

    if (myid==0) {

      if (vm.count("SPL")) psp = std::make_shared<PSPspl>(INITIAL);
      else                 psp = std::make_shared<PSPout>(INITIAL);


      std::cout << "Beginning disk partition [time="
		<< psp->CurrentTime() << "] . . . " << flush;
    }

    partition(psp, compname, particles, histo);
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
  // Open frame list
  // ==================================================

  if (vm.count("SPL")) psp = std::make_shared<PSPspl>(INFILE);
  else                 psp = std::make_shared<PSPout>(INFILE);

  if (myid==0) {
    tnow = psp->CurrentTime();
    std::cout << "Beginning disk partition [time=" << tnow << "] . . . ";
  }

  partition(psp, compname, particles, histo);
  if (myid==0) cout << "done" << endl;

  //------------------------------------------------------------ 

  if (myid==0) cout << "Accumulating for basis . . . " << flush;
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
  write_output(ortho, psp->CurrentTime(), histo);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0) cout << "done" << endl;
  
  //------------------------------------------------------------ 

  MPI_Finalize();

  return 0;
}

