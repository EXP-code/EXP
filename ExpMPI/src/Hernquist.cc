#include "expand.h"
#include <Hernquist.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

Hernquist::Hernquist(string& line) : SphericalBasis(line)
{
  initialize();
  setup();
}

void Hernquist::initialize(void)
{
				// Do nothing . . .
}

double Hernquist::knl(int n, int l)
{
  return 0.5*n*(n+4*l+3) + (l+1)*(2*l+1);
}

double Hernquist::norm(int n, int l)
{
  /*
  return M_PI * knl(n,l) * 
    exp( 
	-log(2.0)*((double)(8*l+4))
	- lgamma((double)(1+n)) - 2.0*lgamma((double)(1.5+2.0*l))
	+ lgamma((double)(4*l+n+3))
	)/(double)(2*l+n+1.5);
  */
  extern double dgammln(double);
  return M_PI * knl(n,l) * 
    exp( 
	-log(2.0)*((double)(8*l+4))
	- dgammln((double)(1+n)) - 2.0*dgammln((double)(1.5+2.0*l))
	+ dgammln((double)(4*l+n+3))
	)/(double)(2*l+n+1.5);
}

void Hernquist::get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp,
			  int tid)
{
  double x, dx, fac, rfac, drfac, dfac1, dfac2;
  int l, n;

  x = r_to_rq(r);
  dx = d_r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfac = 0.5*(1.0 - x);
  drfac = -1.0/(1.0 - x*x);
  
  for (l=0; l<=lmax; l++) {
    dfac1 = 1.0 + x + 2.0*x*l;
    dfac2 = 4.0*l + 3.0;

    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    get_ultra(nmax-1, 2.0*l+1.5, x, du[tid]);

    for (n=0; n<nmax; n++) p[l][n+1] = rfac*u[tid][n];
    dp[l][1] = dx*drfac*dfac1*rfac*u[tid][0];
    for (n=1; n<nmax; n++) dp[l][n+1] = dx*rfac*(drfac*dfac1*u[tid][n] + dfac2*du[tid][n-1]);

    rfac *= fac;
  }

}

void Hernquist::get_potl(int lmax, int nmax, double r, Matrix& p, int tid)
{
  double x, fac, rfac;
  int l, n;

  x = r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfac = 0.5*(1.0 - x);
  
  for (l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    for (n=0; n<nmax; n++) p[l][n+1] = rfac*u[tid][n];
    rfac *= fac;
  }

/* check */
/*  {
    double tmp, potl_HERNQ(int, int, double);

    printf("\nx=%e\n",x);
    for (l=0; l<=lmax; l++) {
      for (n=1; n<nmax; n++) {
	tmp = potl_HERNQ(n, l, x);
	printf("%d %d> %e %e\n",l, n, tmp, p[l][n]);
      }
    }
  } */

}


void Hernquist::get_dens(int lmax, int nmax, double r, Matrix& p, int tid)
{ 
  double x, fac, rfac;
  int l, n;

  x = r_to_rq(r);
  fac = 0.25*(1.0 - x*x);
  rfac = 0.25*pow(1.0 - x, 5.0)/(1.0 - x*x);

  for (l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    for (n=0; n<nmax; n++) p[l][n+1] = krnl[l][n+1]*rfac*u[tid][n];
    rfac *= fac;
  }

}

void Hernquist::get_potl_dens(int lmax, int nmax, double r, 
			      Matrix& p, Matrix& d, int tid)
{
  double x, fac, rfacp, rfacd;
  int l, n;

  x = r_to_rq(r);

  fac = 0.25*(1.0 - x*x);
  rfacp = 0.5*(1.0 - x);
  rfacd = 0.25*pow(1.0 - x, 5.0)/(1.0 - x*x);


  for (l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    for (n=0; n<nmax; n++) p[l][n+1] = rfacp*u[tid][n];
    for (n=0; n<nmax; n++) d[l][n+1] = krnl[l][n+1]*rfacd*u[tid][n];
    rfacp *= fac;
    rfacd *= fac;
  }
}

/*------------------------------------------------------------------------
 *                                                                       *
 *      Convert between reduced coordinate                               *
 *                                                                       *
 *                 r - 1                                                 *
 *          rq =  -------                                                *
 *                 r + 1                                                 *
 *                                                                       *
 *      and its inverse:                                                 *
 *                                                                       *
 *              (1+rq)
 *          r = ------
 *              (1-rq)
 *                                                                       *
 *-----------------------------------------------------------------------*/


#define BIG 1.0e30
double Hernquist::rq_to_r(double rq)
{
  if (rq>=1.0) 
    return BIG;
  else
    return (1.0+rq)/(1.0-rq);
}

double Hernquist::d_r_to_rq(double r)
{
  double fac;

  fac = r + 1.0;
  return 2.0/(fac*fac);
}

double Hernquist::r_to_rq(double r)
{
  return (r-1.0)/(r+1.0);
}
#undef BIG
