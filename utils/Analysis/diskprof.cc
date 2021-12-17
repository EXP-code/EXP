/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute grid slices, and compute volume
 *  for rendering
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
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <string>

#include <yaml-cpp/yaml.h>	// YAML support

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// EXP classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <EmpCylSL.H>
#include <global.H>
#include <localmpi.H>
#include <foarray.H>
#include <DataGrid.H>
#include <writePVD.H>
#include <EXPini.H>

#ifdef DEBUG
#ifndef _REDUCED
#pragma message "NOT using reduced particle structure"
#endif
#endif

const std::string overview = "Compute disk potential, force and density profiles from\nphase-space output files";

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

// Center offset
//
std::vector<double> c0 = {0.0, 0.0, 0.0};

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
    N = std::max<int>(N, 2);
    M = std::max<int>(M, 2);
    dR = 2.0*R/(N-1);
    dZ = 2.0*Z/(M-1);

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
    /*
    dataXY[indX*N + indY] += m;
    dataXZ[indX*M + indZ] += m;
    dataYZ[indY*M + indZ] += m;
    dataZ [indX*N + indY].push_back(z);
    */
  }
};

void add_particles(PRptr reader, const std::string name, vector<Particle>& p, Histogram& h)
{
  // Request particle type
  //
  reader->SelectType(name);

  // Begin reading particles
  //
  auto part = reader->firstParticle();
    
  while (part) {

    // Copy the Particle
    //
    p.push_back(*part);
    
    // Add to histogram
    //
    if (part) h.Add(part->pos[0] - c0[0],
		    part->pos[1] - c0[1],
		    part->pos[2] - c0[2],
		    part->mass);

    // Iterate
    //
    part = reader->nextParticle();
  }

  // Synchronize histogram
  //
  h.Syncr();
}


// Find peak density
// -----------------
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
  double del  = vv[ipeak+1] - vv[ipeak-1];
  double ddel = vv[ipeak+1] + vv[ipeak-1] - 2.0*vv[ipeak];

  if (fabs(ddel) < 1.0e-4) 
    delta = 0.0;		// We're at the peak!
  else
    delta = - 0.5*del/ddel;

  return -ZMAX + dz*(ipeak + delta);

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


void write_output(EmpCylSL& ortho, int icnt, double time, Histogram& histo)
{
  unsigned ncnt = 0;
  int noutV = 9, noutS = 12;
  
  // ==================================================
  // Setup for output files
  // ==================================================
  
  ostringstream sstr;
  sstr << "." << std::setfill('0') << std::setw(5) << icnt;

  string suffix[12] = {"p0", "p1", "p", "fr", "fz", "fp", "d0", "d1", "d",
		       "z10", "z50", "z90"};

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

      std::string OUTF = outdir + "/" + runtag + "_" + outid + "_profile" + sstr.str();
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
    Eigen::VectorXd vv(OUTZ);
    Eigen::VectorXd q;
    std::vector<double> indat(3*OUTR*OUTR, 0.0), otdat(3*OUTR*OUTR);
    
    for (int j=0; j<OUTR; j++) {
	
      x = -RMAX + dR*j;
	
      for (int l=0; l<OUTR; l++) {
      
	y = -RMAX + dR*l;
      
	if ( (ncnt++)%numprocs == myid ) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  for (int k=0; k<OUTZ; k++)
	    vv[k] = ortho.accumulated_dens_eval(r, -ZMAX + dz*k, phi, d0);
	  
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

      DataGrid vtk(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

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
      sout << outdir + "/" + runtag + "_" + outid + "_posn" + sstr.str();
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
    double p0, d0, p1, fr, fz, fp;
    
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
	  
	  ortho.accumulated_eval(r, z, phi, p0, p1, fr, fz, fp);
	  v = ortho.accumulated_dens_eval(r, z, phi, d0);
	  
	  size_t indx = (OUTR*j + l)*OUTR + k;

	  indat[0*blSiz + indx] = p0;
	  indat[1*blSiz + indx] = p1;
	  indat[2*blSiz + indx] = p0 + p1;
	  indat[3*blSiz + indx] = fr;
	  indat[4*blSiz + indx] = fz;
	  indat[5*blSiz + indx] = fp;
	  indat[6*blSiz + indx] = d0;
	  indat[7*blSiz + indx] = v;
	  indat[8*blSiz + indx] = d0 + v;
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
      sout << outdir + "/" + runtag +"_" + outid + "_volume" + sstr.str();
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
    double p0, d0, p1, fr, fz, fp;
    
    vector<double> indat(noutS*OUTR*OUTR, 0.0), otdat(noutS*OUTR*OUTR);
    
    for (int j=0; j<OUTR; j++) {
	
      x = -RMAX + dR*j;
      
      for (int l=0; l<OUTR; l++) {
      
	y = -RMAX + dR*l;
      
	if ((ncnt++)%numprocs == myid) {
	  
	  r = sqrt(x*x + y*y);
	  phi = atan2(y, x);
	  
	  ortho.accumulated_eval(r, z, phi, p0, p1, fr, fz, fp);
	  v = ortho.accumulated_dens_eval(r, z, phi, d0);
	  
	  indat[(0*OUTR+j)*OUTR+l] = p0;
	  indat[(1*OUTR+j)*OUTR+l] = p1;
	  indat[(2*OUTR+j)*OUTR+l] = p0 + p1;
	  indat[(3*OUTR+j)*OUTR+l] = fr;
	  indat[(4*OUTR+j)*OUTR+l] = fz;
	  indat[(5*OUTR+j)*OUTR+l] = fp;
	  indat[(6*OUTR+j)*OUTR+l] = d0;
	  indat[(7*OUTR+j)*OUTR+l] = v;
	  indat[(8*OUTR+j)*OUTR+l] = d0 + v;
	}
      }
    }
    
    MPI_Reduce(&indat[0], &otdat[0], noutS*OUTR*OUTR,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (myid==0) {
      
      for (int j=0; j<OUTR; j++) {
	for (int l=0; l<OUTR; l++) {
	  otdat[(9 *OUTR+j)*OUTR+l] = 0.0;
	  otdat[(10*OUTR+j)*OUTR+l] = 0.0;
	  otdat[(11*OUTR+j)*OUTR+l] = 0.0;
	  
	  // Check for number in the histogram bin
	  //
	  int numZ = histo.dataZ[j*OUTR+l].size();
	  if (numZ>0) {
	    otdat[(9 *OUTR+j)*OUTR+l] = histo.dataZ[j*OUTR+l][floor(0.1*numZ)];
	    otdat[(10*OUTR+j)*OUTR+l] = histo.dataZ[j*OUTR+l][floor(0.5*numZ)];
	    otdat[(11*OUTR+j)*OUTR+l] = histo.dataZ[j*OUTR+l][floor(0.9*numZ)];
	  }
	}
      }

      DataGrid vtkXY(OUTR, OUTR, 1, -RMAX, RMAX, -RMAX, RMAX, 0, 0);

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
      
      DataGrid vtkZ(OUTR, OUTZ, 1, -RMAX, RMAX, -ZMAX, ZMAX, 0, 0);

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
    double p0, d0, p1, fr, fz, fp;
    
    std::vector<double> indat(noutV*OUTR*OUTR, 0.0), otdat(noutV*OUTR*OUTR);
    
      for (int j=0; j<OUTR; j++) {
	
	x = -RMAX + dR*j;

	for (int l=0; l<OUTZ; l++) {
      
	  z = -ZMAX + dZ*l;
      
	  if ((ncnt++)%numprocs == myid) {
	  
	  r   = sqrt(x*x + y*y);
	  phi = atan2(y, x);

	  ortho.accumulated_eval(r, z, phi, p0, p1, fr, fz, fp);
	  v = ortho.accumulated_dens_eval(r, z, phi, d0);
	  
	  indat[(0*OUTR+j)*OUTZ+l] = p0;
	  indat[(1*OUTR+j)*OUTZ+l] = p1;
	  indat[(2*OUTR+j)*OUTZ+l] = p0 + p1;
	  indat[(3*OUTR+j)*OUTZ+l] = fr;
	  indat[(4*OUTR+j)*OUTZ+l] = fz;
	  indat[(5*OUTR+j)*OUTZ+l] = fp;
	  indat[(6*OUTR+j)*OUTZ+l] = d0;
	  indat[(7*OUTR+j)*OUTZ+l] = v;
	  indat[(8*OUTR+j)*OUTZ+l] = d0 + v;
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
      sout << outdir + "/" + runtag + "_" + outid + "_vslice" + sstr.str();
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
    
    for (int j=0; j<OUTR; j++) {
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

	for (int j=0; j<OUTR; j++) {
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

int
main(int argc, char **argv)
{
  int nice, numx, numy, lmax, mmax, nmax, norder, m1, m2, n1, n2;
  int initc, partc, beg, end, stride, init, cmapr, cmapz;
  double rcylmin, rcylmax, rscale, vscale, snr, Hexp=4.0;
  bool DENS, PCA, PVD, verbose = false, mask = false, ignore, logl;
  std::string CACHEFILE, COEFFILE, cname, dir("."), fileType, filePrefix, config;

  // ==================================================
  // MPI preliminaries
  // ==================================================

  local_init_mpi(argc, argv);
  
  // ==================================================
  // Parse command line or input parameter file
  // ==================================================
  
  cxxopts::Options options("diskprof", "Compute disk potential, force and density profiles from phase-space output files");

  options.add_options()
    ("h,help", "produce this help message")
    ("e,expert", "Print the help message showing 'expert' parameters")
    ("v,verbose", "verbose output")
    ("b,mask", "blank empty cells")
    ("T,template", "Write template options file with current and all default values",
     cxxopts::value<string>(config))
    ("f,input", "Input parameter config file",
     cxxopts::value<string>(config))
    ("F,filetype", "input file type",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("P,prefix", "prefix for phase-space files",
     cxxopts::value<std::string>(filePrefix)->default_value("OUT"))
    ("r,runtag", "Runtag name for phase-space files",
     cxxopts::value<string>(runtag)->default_value("run1"))
    ("R,RMAX", "maximum radius for output",
     cxxopts::value<double>(RMAX)->default_value("0.1"))
    ("Z,ZMAX", "maximum height for output",
     cxxopts::value<double>(ZMAX)->default_value("0.01"))
    ("rcylmin", "minimum radius for cylindrical basis table",
     cxxopts::value<double>(rcylmin)->default_value("0.001"))
    ("rcylmax", "maximum radius for cylindrical basis table",
     cxxopts::value<double>(rcylmax)->default_value("20.0"))
    ("NUMX", "number of radial table entries",
     cxxopts::value<int>(numx)->default_value("128"))
    ("NUMY", "number of vertical table entries",
     cxxopts::value<int>(numy)->default_value("64"))
    ("rscale", "radial scale length for basis expansion",
     cxxopts::value<double>(rscale)->default_value("0.01"))
    ("vscale", "vertical scale length for basis expansion",
     cxxopts::value<double>(vscale)->default_value("0.001"))
    ("lmax", "maximum harmonic order for spherical expansion",
     cxxopts::value<int>(lmax)->default_value("36"))
    ("nmax", "maximum harmonic order for spherical expansion",
     cxxopts::value<int>(nmax)->default_value("8"))
    ("mmax", "maximum azimuthal harmonic order for cylindrical expansion",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("norder", "maximum radial order for each harmonic subspace",
     cxxopts::value<int>(norder)->default_value("4"))
    ("outr", "number of radial points for output",
     cxxopts::value<int>(OUTR)->default_value("40"))
    ("outz", "number of vertical points for output",
     cxxopts::value<int>(OUTZ)->default_value("40"))
    ("surface", "make equatorial slices",
     cxxopts::value<bool>(SURFACE)->default_value("true"))
    ("vslice", "make vertical slices",
     cxxopts::value<bool>(VSLICE)->default_value("true"))
    ("probe", "make 1d cuts in and perpendicular to the equatorial plane",
     cxxopts::value<bool>(PROBE)->default_value("true"))
    ("volume", "make volume for rendering",
     cxxopts::value<bool>(VOLUME)->default_value("false"))
    ("axihgt", "compute midplane height profiles",
     cxxopts::value<bool>(AXIHGT)->default_value("false"))
    ("height", "compute height profiles",
     cxxopts::value<bool>(VHEIGHT)->default_value("false"))
    ("C,center", "Accumulation center",
     cxxopts::value<std::vector<double> >(c0))
    ("diff", "render the difference between the trimmed and untrimmed basis")
    ("density", "compute density",
     cxxopts::value<bool>(DENS)->default_value("true"))
    ("compname", "train on Component (default=stars)",
     cxxopts::value<std::string>(cname)->default_value("stars"))
    ("init", "fiducial index",
     cxxopts::value<int>(init)->default_value("0"))
    ("beg", "initial index",
     cxxopts::value<int>(beg)->default_value("0"))
    ("end", "final index",
     cxxopts::value<int>(end)->default_value("99999"))
    ("stride", "index stride",
     cxxopts::value<int>(stride)->default_value("1"))
    ("outdir", "Output directory path",
     cxxopts::value<std::string>(outdir)->default_value("."))
    ("outfile", "Filename prefix",
     cxxopts::value<std::string>(outid)->default_value("diskprof2"))
    ("coeffile", "coefficient output file name",
     cxxopts::value<std::string>(COEFFILE))
    ("cmapr", "Radial coordinate mapping type for cylindrical grid (0=none,, 1=sech, 2=power in z",
     cxxopts::value<int>(cmapr)->default_value("1"))
    ;
  
  options.add_options("expert")
    ("M1", "minimum azimuthal order",
     cxxopts::value<int>(m1)->default_value("0"))
    ("M2", "maximum azimuthal order",
     cxxopts::value<int>(m2)->default_value("1000"))
    ("N1", "minimum radial order",
     cxxopts::value<int>(n1)->default_value("0"))
    ("N2", "maximum radial order",
     cxxopts::value<int>(n2)->default_value("1000"))
    ("nice", "number of bins in x direction",
     cxxopts::value<int>(nice)->default_value("0"))
    ("pvd", "Compute PVD file for ParaView",
     cxxopts::value<bool>(PVD)->default_value("false"))
    ("pca", "perform the PCA analysis for the disk",
     cxxopts::value<bool>(PCA)->default_value("false"))
    ("S,snr", "if not negative: do a SNR cut on the PCA basis",
     cxxopts::value<double>(snr)->default_value("-1.0"))
    ("cachefile", "cachefile name",
     cxxopts::value<std::string>(CACHEFILE)->default_value(".eof.cache.file"))
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
    // Do not overwrite existing config file
    //
    if (std::filesystem::exists(config)) {
      if (myid == 0)
	std::cerr << argv[0] << ": config file <" << config
		  << "> exists, will not overwrite" << std::endl;
      MPI_Finalize();
      return 0;
    }

    // Write template file
    //
    if (myid==0) {
      if (vm.count("expert"))
	SaveConfig(vm, options, config, {"", "expert"});
      else
	SaveConfig(vm, options, config);
    }

    MPI_Finalize();
    return 0;
  }

  // Print help message and exit
  //
  if (vm.count("expert")) {
    if (myid==0) {
      std::cout << std::string(60, '-') << std::endl;
      std::cout << options.help({"", "expert"}) << std::endl;
      std::cout << std::string(60, '-') << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  if (vm.count("help")) {
    if (myid==0) {
      std::cout << std::string(60, '-') << std::endl;
      std::cout << options.help() << std::endl;
      std::cout << std::string(60, '-') << std::endl << std::endl;
    }
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

  if (vm.count("mask")) mask = true;

#ifdef DEBUG
  sleep(20);
#endif  
  
  // ==================================================
  // Nice process
  // ==================================================

  if (nice>0) setpriority(PRIO_PROCESS, 0, nice);

  // ==================================================
  // Parse center
  // ==================================================

  if (vm.count("center")) {
    if (c0.size() != 3) {
      if (myid==0) std::cout << "Center vector needs three components"
			     << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    if (myid==0) {
      std::cout << "Using center: ";
      for (auto v : c0) std::cout << " " << v << " ";
      std::cout << std::endl;
    }
  }

  // ==================================================
  // Phase-space input stream
  // ==================================================

  int iok = 1;

  auto file0 = ParticleReader::fileNameCreator
    (fileType, init, myid, dir, runtag, filePrefix);

  if (myid==0) {
    std::ifstream in(file0);
    if (!in) {
      cerr << "Error opening <" << file0 << ">" << endl;
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

      // Attempt to read magic number
      //
      unsigned int tmagic;
      in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

      //! Basis magic number
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
	DENS    = node["dens"  ].as<bool>();
	if (node["cmap"])
	  cmapr = node["cmap"  ].as<int>();
	else
	  cmapr = node["cmapr" ].as<int>();
	if (node["cmapz"])
	  cmapz = node["cmapz"  ].as<int>();
	rcylmin = node["rmin"  ].as<double>();
	rcylmax = node["rmax"  ].as<double>();
	rscale  = node["ascl"  ].as<double>();
	vscale  = node["hscl"  ].as<double>();
	
      } else {
				// Rewind file
	in.clear();
	in.seekg(0);

	int tmp;
    
	in.read((char *)&mmax,    sizeof(int));
	in.read((char *)&numx,    sizeof(int));
	in.read((char *)&numy,    sizeof(int));
	in.read((char *)&nmax,    sizeof(int));
	in.read((char *)&norder,  sizeof(int));
      
	in.read((char *)&tmp,     sizeof(int)); 
	if (tmp) DENS = true;
	else     DENS = false;
	
	in.read((char *)&cmapr,   sizeof(int)); 
	in.read((char *)&rcylmin, sizeof(double));
	in.read((char *)&rcylmax, sizeof(double));
	in.read((char *)&rscale,  sizeof(double));
	in.read((char *)&vscale,  sizeof(double));
      }
    }
  }

  EmpCylSL::RMIN        = rcylmin;
  EmpCylSL::RMAX        = rcylmax;
  EmpCylSL::NUMX        = numx;
  EmpCylSL::NUMY        = numy;
  EmpCylSL::CMAPR       = cmapr;
  EmpCylSL::CMAPZ       = cmapz;
  EmpCylSL::logarithmic = logl;
  EmpCylSL::DENS        = DENS;
  EmpCylSL::CACHEFILE   = CACHEFILE;

				// Create expansion
				//
  EmpCylSL ortho(nmax, lmax, mmax, norder, rscale, vscale);

				// Range limits for evaluation
				//
  ortho.set_mrange(m1, m2);
  ortho.set_nrange(n1, n2);

  Histogram histo(OUTR, OUTZ, RMAX, ZMAX);
  std::vector<Particle> particles;
  PRptr reader;

  std::vector<double> times;
  std::vector<std::string> outfiles;

  bool compute = false;
  if (PCA) {
    EmpCylSL::PCAVAR = true;
    EmpCylSL::HEXP   = Hexp;
    if (vm.count("truncate"))
      ortho.setTK("Truncate");
    else
      ortho.setTK("Hall");
    
    compute = true;
    ortho.setup_accumulation();
  }

  // ==================================================
  // Initialize and/or create basis
  // ==================================================
  
  if (ortho.read_cache()==0) {
    
    try {
      reader = ParticleReader::createReader(fileType, file0, myid, true);
    }
    catch (std::runtime_error &error) {
      std::cerr << error.what() << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    if (myid==0) std::cout << "Beginning disk partition [time="
			   << reader->CurrentTime()
			   << "] . . . " << std::flush;
      
    particles.clear();		// Empty the particle vector
				// Then read in new particles
    add_particles(reader, cname, particles, histo);

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
    ortho.make_coefficients(compute);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
  }
  

  // ==================================================
  // Open output coefficient file
  // ==================================================
  
  std::ofstream coefs;
  if (myid==0 and COEFFILE.size()>0) {
    COEFFILE = outdir + "/" + COEFFILE;
    coefs.open(COEFFILE);
    if (not coefs) {
      std::cerr << "Could not open coefficient file <" << COEFFILE << "> . . . quitting"
		<< std::endl;
      MPI_Abort(MPI_COMM_WORLD, 10);
    }
  }

  std::string file1;

  for (int indx=beg; indx<=end; indx+=stride) {

    // ==================================================
    // Phase-space reader input stream
    // ==================================================

    iok = 1;

    file1 = ParticleReader::fileNameCreator
      (fileType, indx, myid, dir, runtag, filePrefix);

    if (myid==0) {
      std::ifstream in(file1);
      if (!in) {
	cerr << "Error opening <" << file1 << ">" << endl;
	iok = 0;
      }
    }
    
    MPI_Bcast(&iok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (iok==0) break;
    
    // ==================================================
    // Open frame list
    // ==================================================
    
    try {
      reader = ParticleReader::createReader(fileType, file1, myid, true);
    }
    catch (std::runtime_error &error) {
      std::cerr << error.what() << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    double tnow = reader->CurrentTime();

    if (myid==0) {
      cout << "Beginning disk partition [time=" << tnow
	   << ", index=" << indx << "] . . . "  << flush;

      times.push_back(tnow);
      ostringstream filen;
      filen << runtag << "_" << outid << "_surface."
	    << std::setfill('0') << std::setw(5) << indx << ".vtr";
      outfiles.push_back(filen.str());
    }
      
    histo.Reset();		// Reset surface histogram

    particles.clear();		// Empty the particle vector
				// Then read in new particles
    add_particles(reader, cname, particles, histo);
    if (myid==0) cout << "done" << endl;
    
    ortho.setup_accumulation();
    if (PCA)
      ortho.setHall(outdir + "/" + runtag + "_" + outid + ".pca", particles.size());

    if (myid==0) cout << "Accumulating particle positions . . . " << flush;
    ortho.accumulate(particles, 0, true, compute);
    //                             ^     ^
    //                             |     |
    // Verbose --------------------+     |
    // Compute covariance ---------------+
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
      
    //------------------------------------------------------------ 
      
    if (myid==0) cout << "Making disk coefficients . . . " << flush;
    ortho.make_coefficients(compute);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
      
    //------------------------------------------------------------ 

    if (myid==0) {
      std::cout << "Writing disk coefficients . . . " << flush;
      ortho.dump_coefs_binary(coefs, tnow);
      cout << "done" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //------------------------------------------------------------ 
      
    if (myid==0) cout << "Writing output . . . " << flush;
    //
    // Coefficient trimming
    //
    if (PCA and snr>=0.0) {
      std::vector<Eigen::VectorXd> ac_cos, ac_sin;
      std::vector<Eigen::VectorXd> rt_cos, rt_sin, sn_rat;
      ortho.pca_hall(true, false);
      ortho.get_trimmed(snr, ac_cos, ac_sin,
			&rt_cos, &rt_sin, &sn_rat);
      for (int mm=0; mm<=mmax; mm++) {
	if (vm.count("diff")) {
	  ac_cos[mm] = ac_cos[mm] - rt_cos[mm];
	  if (mm) 
	    ac_sin[mm] = ac_sin[mm] - rt_sin[mm];
	}
	if (mm==0)
	  ortho.set_coefs(mm, ac_cos[mm], ac_sin[mm], true);
	else
	  ortho.set_coefs(mm, ac_cos[mm], ac_sin[mm], false);
      }
      // END: M loop
    }

    double time = 0.0;
    if (myid==0) time = reader->CurrentTime();
    write_output(ortho, indx, time, histo);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done" << endl;
    
    //------------------------------------------------------------ 
    
  } // Dump loop

#ifdef HAVE_VTK
  // Create PVD file
  //
  if (myid==0 and PVD) {
    writePVD(outdir + "/" + runtag + ".pvd", times, outfiles);
  }
#endif

  // Shutdown MPI
  //
  MPI_Finalize();
    
  // DONE
  //
  return 0;
  }

