// This is really C++ code
/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the potential, acceleration and density using
 *  periodic box expansion in X & Y and slab in Z
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
#include <biorth1d.h>

extern "C" {
#include "expand.h"
}

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

static const double KEPS=1.0e-6;
static Complex *expccof;
static double *expreal, *expimag;
static double *expreal1, *expimag1;
static OneDTrig **trig;

struct SlabCoefHeader {
  double time;
  double zmax;
  int nmaxx, nmaxy, nmaxz;
  int jmax;
};

static SlabCoefHeader coefheader;

void 
acceleration_and_potential_slab(void)
{
  static int firstime=1;
  static int imx,imy,imz,jmax,nnmax;
  int i, j, ix, iy, iz, iix, iiy, ii, jj, indx, dum;
  static double  dfac;
  static Complex kfac;
  Complex fac, startx, starty, facx, facy, potl, facf;
  Complex accx, accy, accz;
  Complex stepx, stepy;
  static Vector zfrc, zpot;

  if (firstime) {
    firstime=0;

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

    trig = new OneDTrig* [nnmax+1];
    for (i=0; i<=nnmax; i++) {
      trig[i] = new OneDTrig [i+1];
      for (j=0; j<=i; j++) trig[i][j].reset(dfac*sqrt(i*i + j*j)+KEPS, zmax);
    }

    zpot.setsize(1, nmaxz);
    zfrc.setsize(1, nmaxz);
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
    expreal1[indx] = 0.0;
    expimag1[indx] = 0.0;
  }

  
  /*            Compute coefficients */

  if (myid) {

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
	    trig[iix][iiy].potl(dum, dum, z[i+1], zpot);
	  else
	    trig[iiy][iix].potl(dum, dum, z[i+1], zpot);

	  for (iz=0; iz<imz; iz++) {

	    indx = imz*(iy + imy*ix) + iz;
	    expccof[indx] += 4.0*M_PI*mass[i+1]*facx*facy*zpot[iz+1];
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

  }

  MPI_Allreduce( expreal1, expreal, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce( expimag1, expimag, jmax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (indx=0; indx<jmax; indx++)
    expccof[indx] = Complex(expreal[indx], expimag[indx]);


  /* Determine potential and acceleration */

  if (myid) {

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
	    trig[iix][iiy].potl(dum, dum, z[i+1], zpot);
	    trig[iix][iiy].force(dum, dum, z[i+1], zfrc);
	  }
	  else {
	    trig[iiy][iix].potl(dum, dum, z[i+1], zpot);
	    trig[iiy][iix].force(dum, dum, z[i+1], zfrc);
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

	    potl -= fac;
	    
	    accx -= Complex(0.0,-dfac*ii)*fac;
	    accy -= Complex(0.0,-dfac*jj)*fac;
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

}

#include <stdio.h>

extern "C" void dump_coefs_slab(FILE *fout)
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


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------

void usage(char *);
void set_parm(char *, char *);
void write_parm(void);
void print_parm(ostream &, char *);
void print_default(void);

extern "C" {

  void get_acceleration_and_potential_slab(void)
    {
      acceleration_and_potential_slab();
    }
}

