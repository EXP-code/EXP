/*****************************************************************************
 *  Description:
 *  -----------
 *
 *
 *  Routines for computing biorthonormal pairs for the Hernquist model 
 *  based on Clutton-Brock's ultraspherical series
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

double r_to_rq(double r),rq_to_r(double rq),normHERNQ(int n, int l),knlHERNQ(int n, int l);

double potl_HERNQ(int nn, int l, double x)
{
  double ultra(int n, double dl, double x);
  int n;

/*  x = r_to_rq(r); */

  if (fabs(x) >= 1.0) return 0.0;

  n = nn-1;
  return pow(1.0 - x*x,(double)l) * (1.0 - x)* ultra(n,2.0*l+0.5,x)/
    pow(2.0,2.0*l + 1.0);

}



double dens_HERNQ(int nn, int l, double x)
{
  double ultra(int n, double dl, double x),dgammln();
  int n;

/*  x = r_to_rq(r); */

  if (fabs(x) >= 1.0) return 0.0;

  n = nn-1;
  return knlHERNQ(n,l)/pow(2.0,2.0*l + 2.0) *
    pow(1.0 - x*x,(double)l-1.0) * pow(1.0 - x,5.0) * ultra(n,2.0*l+0.5,x);
}

double knlHERNQ(int n, int l)
{
  return 0.5*n*(n+4*l+3) + (l+1)*(2*l+1);
}

double normHERNQ(int n, int l)
{
  return M_PI * knlHERNQ(n,l) * exp( 
	    -log(2.0)*((double)(8*l+4))
	    - dgammln((double)(1+n)) - 2.0*dgammln((double)(1.5+2.0*l))
		    + dgammln((double)(4*l+n+3))
		    )/(double)(2*l+n+1.5);
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
double rq_to_r(double rq)
{
  if (rq>=1.0) 
    return BIG;
  else
    return (1.0+rq)/(1.0-rq);
}

double d_r_to_rq(double r)
{
  double fac;

  fac = r + 1.0;
  return 2.0/(fac*fac);
}

double r_to_rq(double r)
{
  return (r-1.0)/(r+1.0);
}
#undef BIG


