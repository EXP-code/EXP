// This may look like C code, but it is really -*- C++ -*-

// #define DEBUG 1

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  the SL biorthogonal expansion
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
 *  Value
 *
 *  Notes:
 *  -----
 *
 *  By:
 *  --
 *
 *  MDW 03/24/98
 *
 ***************************************************************************/

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

#define NGRID 2000

#include <stdlib.h>
#include <math.h>
#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <kevin_complex.h>
#include <SLGridMP2.h>

void determine_coefficients_slabSL(void);
void determine_acceleration_and_potential_slabSL(void);

//======================================================================
//======================================================================

static const double KEPS=1.0e-6;
static Complex *expccof;
static double *expreal, *expimag;
static double *expreal1, *expimag1;
static int imx, imy, imz, jmax;
static double  dfac;
static Complex kfac;
static Vector zfrc, zpot;

static SLGridSlab *grid;

struct SlabCoefHeader {
  double time;
  double zmax;
  int nmaxx, nmaxy, nmaxz;
  int jmax;
};
static SlabCoefHeader coefheader;

//======================================================================
//======================================================================

extern "C" 
void get_acceleration_and_potential_slabSL(void)
{
  static int firstime=1;

  if (firstime) {

    SLGridSlab::mpi = 1;
    SLGridSlab::ZBEG = 0.0;
    SLGridSlab::ZEND = 0.1;
    SLGridSlab::H = hslab;
  
    int nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

    grid = new SLGridSlab(nnmax, nmaxz, NGRID, zmax);

    imx = 1+2*nmaxx;
    imy = 1+2*nmaxy;
    imz = nmaxz;
    jmax = imx*imy*imz;
    expccof = new Complex[jmax];

    expreal = new double [jmax];
    expreal1 = new double [jmax];
    expimag = new double [jmax];
    expimag1 = new double [jmax];
    
    dfac = 2.0*M_PI;
    kfac = Complex(0.0, dfac);

    nnmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;

    zpot.setsize(1, nmaxz);
    zfrc.setsize(1, nmaxz);
  }

  /*======================*/
  /* Compute coefficients */
  /*======================*/

  if (firstime || self_consistent) {
    firstime = 0;
    determine_coefficients_slabSL();
  }

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

  determine_acceleration_and_potential_slabSL();

  MPL_stop_timer();

}


void 
determine_coefficients_slabSL(void)
{
  int i, ix, iy, iz, iix, iiy, indx;
  Complex startx, starty, facx, facy;
  Complex stepx, stepy;



  /*======================*/
  /* Compute coefficients */
  /*======================*/

				//  Coefficients are ordered as follows:
				//  n=-nmax,-nmax+1,...,0,...,nmax-1,nmax
				//  in a single array for each dimension
				//  with z dimension changing most rapidly
  /*		Clean */
  for (indx=0; indx<jmax; indx++) {
    expccof[indx] = 0.0;
    expreal[indx] = 0.0;
    expimag[indx] = 0.0;
    expreal1[indx] = 0.0;
    expimag1[indx] = 0.0;
  }

  
  /*            Compute coefficients */


  for (i=0; i<nbodies; i++) {

				/* Truncate to box with sides in [0,1] */
    
    if (x[i+1]<0.0)
      x[i+1] += (double)((int)fabs(x[i+1])) + 1.0;
    else
      x[i+1] -= (double)((int)x[i+1]);
    
    if (y[i+1]<0.0)
      y[i+1] += (double)((int)fabs(y[i+1])) + 1.0;
    else
      y[i+1] -= (double)((int)y[i+1]);
    
    
				/* Recursion multipliers */
    stepx = exp(-kfac*x[i+1]);
    stepy = exp(-kfac*y[i+1]);
    
				/* Initial values */
    startx = exp(nmaxx*kfac*x[i+1]);
    starty = exp(nmaxy*kfac*y[i+1]);
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
      iix = nmaxx - ix;
      if (iix<0) iix = ix - nmaxx;
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	iiy = nmaxy - iy;
	if (iiy<0) iiy = iy - nmaxy;
	
	if (iix < 0 || iix > nmaxx) {
	  cerr << "Out of bounds: iix=" << iix << endl;
	}
	if (iiy < 0 || iiy > nmaxy) {
	  cerr << "Out of bounds: iiy=" << iiy << endl;
	}
	
	if (iix>=iiy)
	  grid->get_pot(zpot, z[i+1], iix, iiy);
	else
	  grid->get_pot(zpot, z[i+1], iiy, iix);

	for (iz=0; iz<imz; iz++) {

	  indx = imz*(iy + imy*ix) + iz;
	  expccof[indx] += -4.0*M_PI*mass[i+1]*facx*facy*zpot[iz+1];
				// ^
				// |--- density in orthogonal series
				//      is 4.0*M_PI rho
	}
      }
    }
  }
    
  
  for (indx=0; indx<jmax; indx++) {
    expreal1[indx] = expccof[indx].real();
    expimag1[indx] = expccof[indx].imag();
  }
  

  MPI_Allreduce( expreal1, expreal, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce( expimag1, expimag, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (indx=0; indx<jmax; indx++)
    expccof[indx] = Complex(expreal[indx], expimag[indx]);

}


void
determine_acceleration_and_potential_slabSL(void)
{
  int i, ix, iy, iz, iix, iiy, ii, jj, indx;
  Complex fac, startx, starty, facx, facy, potl, facf;
  Complex stepx, stepy;
  Complex accx, accy, accz;


  /* Determine potential and acceleration */

  
  for (i=0; i<nbodies; i++) {
    
    accx = accy = accz = potl = 0.0;
    
				/* Recursion multipliers */
    stepx = exp(kfac*x[i+1]);
    stepy = exp(kfac*y[i+1]);
    
				/* Initial values (note sign change) */
    startx = exp(-nmaxx*kfac*x[i+1]);
    starty = exp(-nmaxy*kfac*y[i+1]);
    
    for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
      
				/* Compute wavenumber; recall that the */
				/* coefficients are stored as follows: */
				/* -nmax,-nmax+1,...,0,...,nmax-1,nmax */
      ii = ix - nmaxx;
      iix = abs(ii);
      
      for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	
	jj = iy - nmaxy;
	iiy = abs(jj);
	
	if (iix < 0 || iix > nmaxx) {
	  cerr << "Out of bounds: iix=" << iix << endl;
	}
	if (iiy < 0 || iiy > nmaxy) {
	  cerr << "Out of bounds: iiy=" << iiy << endl;
	}
	
	if (iix>=iiy) {
	  grid->get_pot  (zpot, z[i+1], iix, iiy);
	  grid->get_force(zfrc, z[i+1], iix, iiy);
	}
	else {
	  grid->get_pot  (zpot, z[i+1], iiy, iix);
	  grid->get_force(zfrc, z[i+1], iiy, iix);
	}
	
	for (iz=0; iz<imz; iz++) {
	  
	  indx = imz*(iy + imy*ix) + iz;
	  
	  fac  = facx*facy*zpot[iz+1]*expccof[indx];
	  facf = facx*facy*zfrc[iz+1]*expccof[indx];
	  
	  
				/* Don't add potential for 
				   zero wavenumber (constant term) */
	  
	  if (ii==0 && jj==0 && iz==0) fac = 0.0;
	  
	  
				/* Limit to minimum wave number */
	  
	  if (abs(ii)<nminx || abs(jj)<nminy) continue;
	  
	  potl += fac;
	  
	  accx += Complex(0.0,-dfac*ii)*fac;
	  accy += Complex(0.0,-dfac*jj)*fac;
	  accz -= facf;
	  
	}
      }
    }
    
    ax[i+1] += Re(accx);
    ay[i+1] += Re(accy);
    az[i+1] += Re(accz);
    pot[i+1] += Re(potl);
  }
  
}

#include <stdio.h>

extern "C" void dump_coefs_slabSL(FILE *fout)
{
  coefheader.time = tnow;
  coefheader.zmax = zmax;
  coefheader.nmaxx = nmaxx;
  coefheader.nmaxy = nmaxy;
  coefheader.nmaxz = nmaxz;
  coefheader.jmax = (1+2*nmaxx)*(1+2*nmaxy)*nmaxz;
  
  fwrite(&coefheader, sizeof(SlabCoefHeader), 1, fout);
  for (int i=0; i<coefheader.jmax; i++) {
    fwrite(&expccof[i].real(), sizeof(double), 1, fout);
    fwrite(&expccof[i].imag(), sizeof(double), 1, fout);
  }
}



