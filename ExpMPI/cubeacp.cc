// This is really C++ code
/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  periodic cube expansion
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *  G=1
 *
 *  By:
 *  --
 *
 *  MDW 04/27/92; revision 04/29/92
 *
 *  The MPI revision is untested
 *
 ***************************************************************************/
#define CPLUSPLUS

#include <stdlib.h>
#include <math.h>
#include <kevin_complex.h>

extern "C" {
#include "expand.h"
}

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

static Complex *expccof;
static double *expreal, *expimag;
static double *expreal1, *expimag1;


void acceleration_and_potential_cube(void)
{
  static int firstime=1;
  static int imx,imy,imz,jmax;
  int i,ix,iy,iz,ii,jj,kk,indx;
  static double  dfac;
  static Complex kfac;
  Complex fac,startx,starty,startz,facx,facy,facz,dens,potl,accx,accy,accz;
  Complex stepx,stepy,stepz;
  double k2;

  if (firstime) {
    firstime=0;

    imx = 1+2*nmaxx;
    imy = 1+2*nmaxy;
    imz = 1+2*nmaxz;
    jmax = imx*imy*imz;
    expccof = new Complex[jmax];

    expreal = new double [jmax];
    expreal1 = new double [jmax];
    expimag = new double [jmax];
    expimag1 = new double [jmax];
    
    dfac = 2.0*M_PI;
    kfac = Complex(0.0,dfac);
  }


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
  }

  
  /*            Compute coefficients */

  if (myid) {

    for (i=0; i<nbodies; i++) {

				/* Truncate to cube with sides in [0,1] */

      if (x[i+1]<0.0)
	x[i+1] += (double)((int)fabs(x[i+1])) + 1.0;
      else
	x[i+1] -= (double)((int)x[i+1]);

      if (y[i+1]<0.0)
	y[i+1] += (double)((int)fabs(y[i+1])) + 1.0;
      else
	y[i+1] -= (double)((int)y[i+1]);

      if (z[i+1]<0.0)
	z[i+1] += (double)((int)fabs(z[i+1])) + 1.0;
      else
	z[i+1] -= (double)((int)z[i+1]);


				/* Recursion multipliers */
      stepx = exp(-kfac*x[i+1]);
      stepy = exp(-kfac*y[i+1]);
      stepz = exp(-kfac*z[i+1]);

				/* Initial values */
      startx = exp(nmaxx*kfac*x[i+1]);
      starty = exp(nmaxy*kfac*y[i+1]);
      startz = exp(nmaxz*kfac*z[i+1]);

      for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
	for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	  for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {

	    indx = imz*(iy + imy*ix) + iz;
	    expccof[indx] += mass[i+1]*facx*facy*facz;
	  
	  }
	}
      }
    }
    

    for (indx=0; indx<jmax; indx++) {
      expreal1[indx] = expccof[indx].real();
      expimag1[indx] = expccof[indx].imag();
    }

  }

  MPI_Allreduce( expreal1, expreal, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce( expimag1, expimag, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);

  for (indx=0; indx<jmax; indx++)
    expccof[indx] = Complex(expreal[indx], expimag[indx]);


  /* Determine potential and acceleration */

  if (myid) {

    for (i=0; i<nbodies; i++) {

      accx = accy = accz = dens = potl = 0.0;

				/* Recursion multipliers */
      stepx = exp(kfac*x[i+1]);
      stepy = exp(kfac*y[i+1]);
      stepz = exp(kfac*z[i+1]);

				/* Initial values (note sign change) */
      startx = exp(-nmaxx*kfac*x[i+1]);
      starty = exp(-nmaxy*kfac*y[i+1]);
      startz = exp(-nmaxz*kfac*z[i+1]);

      for (facx=startx, ix=0; ix<imx; ix++, facx*=stepx) {
	for (facy=starty, iy=0; iy<imy; iy++, facy*=stepy) {
	  for (facz=startz, iz=0; iz<imz; iz++, facz*=stepz) {

	    indx = imz*(iy + imy*ix) + iz;

	    fac = facx*facy*facz*expccof[indx];
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

      ax[i+1] += Re(accx);
      ay[i+1] += Re(accy);
      az[i+1] += Re(accz);
      pot[i+1] += Re(potl);
    }

  }

}

			// C interface

extern "C" {

  void get_acceleration_and_potential_cube(void)
    {
      acceleration_and_potential_cube();
    }
}

