/*****************************************************************************
 *  Description:
 *  -----------
 *
 *
 *  Routines for computing biorthonormal pairs based on 
 *  Clutton-Brock's ultraspherical series
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
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91 [based on "findwake" by MDW: 12/26/1987]
 *
 ***************************************************************************/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

double r_to_rh(double r),rh_to_r(double rh),normCB(int n, int l),knlCB(int n, int l);

double potl_CB(int nn, int l, double x)
{
  double ultra(int n, double dl, double x);
  int n;

/*  x = r_to_rh(r); */

  if (fabs(x) >= 1.0) return 0.0;

  n = nn-1;
  return pow(1.0 - x*x,0.5*(double)l) * sqrt(1.0 - x) * ultra(n,(double)l,x)/
    pow(2.0,0.5+(double)l);
}



double dens_CB(int nn, int l, double x)
{
  double ultra(int n, double dl, double x),dgammln();
  int n;

/*  x = r_to_rh(r); */

  if (fabs(x) >= 1.0) return 0.0;

  n = nn-1;
  return knlCB(n,l)/pow(2.0,2.5+(double)l) *
    pow(1.0 - x*x,0.5*(double)l) * pow(1.0 - x,2.5) * ultra(n,(double)l,x);
}

double knlCB(int n, int l)
{
  return 4.0*n*(n+2*l+2) + (2*l+1)*(2*l+3);
}

double normCB(int n, int l)
{
  return M_PI * knlCB(n,l) * exp( 
		    -log(2.0)*((double)(4*l+4))
		    - dgammln((double)(1+n)) - 2.0*dgammln((double)(1+l) )
		    + dgammln((double)(2*l+n+2))
		    )/(double)(l+n+1);
}


/*------------------------------------------------------------------------
 *                                                                       *
 *      Convert between reduced coordinate                               *
 *                                                                       *
 *                 r^2-1                                                 *
 *          rh =  -------                                                *
 *                 r^2+1                                                 *
 *                                                                       *
 *      and its inverse:                                                 *
 *                                                                       *
 *              (1+rh)^(1/2)                                             *
 *          r = ------------                                             *
 *              (1-rh)^(1/2)                                             *
 *                                                                       *
 *-----------------------------------------------------------------------*/


#define BIG 1.0e30
double rh_to_r(double rh)
{
  if (rh>=1.0) 
    return BIG;
  else
    return sqrt( (1.0+rh)/(1.0-rh) );
}

double d_r_to_rh(double r)
{
  double fac;

  fac = r*r + 1.0;;
  return 4.0*r/(fac*fac);
}

double r_to_rh(double r)
{
  return (r*r-1.0)/(r*r+1.0);
}
#undef BIG

