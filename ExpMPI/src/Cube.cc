#include "expand.h"

#include <stdlib.h>
#include <math.h>

#include <Cube.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

Cube::Cube(string& line) : PotAccel(line)
{
  id = "Cube";

  nminx = nminy = nminz = 0;
  nmaxx = nmaxy = nmaxz = 16;

  initialize();

  imx = 1+2*nmaxx;
  imy = 1+2*nmaxy;
  imz = 1+2*nmaxz;
  jmax = imx*imy*imz;
  expccof = new Complex* [nthrds];
  for (int i=0; i<nthrds; i++)
    expccof[i] = new Complex[jmax];

  expreal = new double [jmax];
  expreal1 = new double [jmax];
  expimag = new double [jmax];
  expimag1 = new double [jmax];
  
  dfac = 2.0*M_PI;
  kfac = Complex(0.0,dfac);
}

Cube::~Cube(void)
{
  for (int i=0; i<nthrds; i++) delete [] expccof[i];
  delete [] expccof;

  delete [] expreal;
  delete [] expreal1;
  delete [] expimag;
  delete [] expimag1;
}

void Cube::initialize(void)
{
  string val;

  if (get_value("nminx", val)) nminx = atoi(val.c_str());
  if (get_value("nminy", val)) nminy = atoi(val.c_str());
  if (get_value("nminz", val)) nminz = atoi(val.c_str());
  if (get_value("nmaxx", val)) nmaxx = atoi(val.c_str());
  if (get_value("nmaxy", val)) nmaxy = atoi(val.c_str());
  if (get_value("nmaxz", val)) nmaxz = atoi(val.c_str());
}

void * Cube::determine_coefficients_thread(void * arg)
{

  int i,ix,iy,iz,indx;
  Complex startx,starty,startz,facx,facy,facz;
  Complex stepx,stepy,stepz;
  double mass;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  use[id] = 0;

  for (i=nbeg; i<nend; i++) {

    use[id]++;
    mass = (*particles)[i].mass;

				/* Truncate to cube with sides in [0,1] */
    
    if ((*particles)[i].pos[0]<0.0)
      (*particles)[i].pos[0] += 
	(double)((int)fabs((*particles)[i].pos[0])) + 1.0;
    else
      (*particles)[i].pos[0] -= (double)((int)(*particles)[i].pos[0]);
    
    if ((*particles)[i].pos[1]<0.0)
      (*particles)[i].pos[1] += 
	(double)((int)fabs((*particles)[i].pos[1])) + 1.0;
    else
      (*particles)[i].pos[1] -= (double)((int)(*particles)[i].pos[1]);
    
    if ((*particles)[i].pos[2]<0.0)
      (*particles)[i].pos[2] += 
	(double)((int)fabs((*particles)[i].pos[2])) + 1.0;
    else
      (*particles)[i].pos[2] -= (double)((int)(*particles)[i].pos[2]);
    
    
				/* Recursion multipliers */
    stepx = exp(-kfac*(*particles)[i].pos[0]);
    stepy = exp(-kfac*(*particles)[i].pos[1]);
    stepz = exp(-kfac*(*particles)[i].pos[2]);
    
				/* Initial values */
    startx = exp(nmaxx*kfac*(*particles)[i].pos[0]);
    starty = exp(nmaxy*kfac*(*particles)[i].pos[1]);
    startz = exp(nmaxz*kfac*(*particles)[i].pos[2]);
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  expccof[id][indx] += mass*facx*facy*facz;
	  
	}
      }
    }
  }
    
  return (NULL);
}

void Cube::determine_coefficients(void)
{
				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly
  // Clean 
  for (int indx=0; indx<jmax; indx++) {
    for (int i=0; i<nthrds; i++) expccof[i][indx] = 0.0;
    expreal[indx] = 0.0;
    expimag[indx] = 0.0;
  }

  exp_thread_fork(true);

  for (int i=0; i<nthrds; i++) use1 += use[i];
  
  MPI_Allreduce ( &use1, &use0,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  used = use0;

  for (int i=0; i<nthrds; i++) {
    for (int indx=0; indx<jmax; indx++) {
      expreal1[indx] = expccof[i][indx].real();
      expimag1[indx] = expccof[i][indx].imag();
    }
  }

  MPI_Allreduce( expreal1, expreal, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce( expimag1, expimag, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int indx=0; indx<jmax; indx++)
    expccof[0][indx] = Complex(expreal[indx], expimag[indx]);

}

void * Cube::determine_acceleration_and_potential_thread(void * arg)
{
  int ix,iy,iz,ii,jj,kk,indx;
  Complex fac,startx,starty,startz,facx,facy,facz,dens,potl,accx,accy,accz;
  Complex stepx,stepy,stepz;
  double k2;

  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (int i=nbeg; i<nend; i++) {
    
    accx = accy = accz = dens = potl = 0.0;
    
				/* Recursion multipliers */
    stepx = exp(kfac*(*particles)[i].pos[0]);
    stepy = exp(kfac*(*particles)[i].pos[1]);
    stepz = exp(kfac*(*particles)[i].pos[2]);
    
				/* Initial values (note sign change) */
    startx = exp(-nmaxx*kfac*(*particles)[i].pos[0]);
    starty = exp(-nmaxy*kfac*(*particles)[i].pos[1]);
    startz = exp(-nmaxz*kfac*(*particles)[i].pos[2]);
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  
	  fac = facx*facy*facz*expccof[0][indx];
	  dens += fac;
	  
				/* Compute wavenumber; recall that the */
				/* coefficients are stored as follows: */
				/* -nmax,-nmax+1,...,0,...,nmax-1,nmax */
	  ii = ix-nmaxx;
	  jj = iy-nmaxy;
	  kk = iz-nmaxz;
	  
				/* No contribution to acceleration and */
	                        /* potential ("swindle") for zero      */
	                        /* wavenumber */
	  
	  if (ii==0 && jj==0 && kk==0) continue;
	  
				/* Limit to minimum wave number */
	  
	  if (abs(ii)<nminx || abs(jj)<nminy || abs(kk)<nminz) continue;
	  
	  k2 = 4.0*M_PI/(dfac*dfac*(ii*ii + jj*jj + kk*kk));
	  potl -= k2*fac;
	  
	  accx -= k2*Complex(0.0,-dfac*ii)*fac;
	  accy -= k2*Complex(0.0,-dfac*jj)*fac;
	  accz -= k2*Complex(0.0,-dfac*kk)*fac;
	  
	}
      }
    }
    
    (*particles)[i].acc[0] += Re(accx);
    (*particles)[i].acc[1] += Re(accy);
    (*particles)[i].acc[2] += Re(accz);
    if (use_external)
      (*particles)[i].potext += Re(potl);
    else
      (*particles)[i].pot += Re(potl);
  }
  
  return (NULL);
}

void Cube::get_acceleration_and_potential(vector<Particle>* Particles)
{
  particles = Particles;

  exp_thread_fork(false);

  // Clear external potential flag
  use_external = false;
}
