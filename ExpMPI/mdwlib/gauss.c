/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine performs an intgrals of the form
 *
 *   Type 0:
 *
 *           b
 *           /
 *       I = | dx f(x)
 *           /
 *           a
 *
 *
 *   Type 2:
 *
 *           b
 *           /     f(x)
 *       I = | dx ------
 *           /         1/2
 *           a    [x-a]
 *
 *  using Gaussian quadrature.
 *
 *  Call sequence:
 *  -------------
 *  x = gint_0(a,b,f,NGauss);
 *  x = gint_2(a,b,f,NGauss);
 *
 *  double x,a,b,f(),gint();
 *  int NGauss;
 *
 *  Parameters:
 *  ----------
 *
 *  a        as above
 *  b        as above
 *  f        function f(x)
 *  NGauss   number of "knots" in quadrature.  If NGauss is changed from
 *           last call, weights and abscissas are recomputed otherwise
 *           set of values are used.
 *
 *  Returns:
 *  -------
 *
 *  Value of integral
 *
 *  Notes:
 *  -----
 *
 *  Routine must be linked with 'cutil' or 'nr' as it uses vector allocations.
 *
 *  11/13/88 Checked.  Values returned by gauss_leg_knots() are agree with
 *          tables in A&S to (at least) 15 places 
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include "cutil.h"
                              /* if DEBUG is set, weights and abscissas are
                                 printed on stderr */
#ifndef DEBUG
#   define DEBUG 0
#endif

static int N=0,M;             /* N=0 force comp. of knots on first call */
static double *x,*w;
void gauss_leg_knots(void);

double gint_0(double a, double b, double (*f) (double), int NGauss)
{
  double accum,bma,bpa;
  int i;

/* Get values weights and abcissas for Gauss-Legendre integration if needed */

  if (NGauss != N) {
    if (N != 0) {
      free_dvector(x,1,M);
      free_dvector(w,1,M);
    }
    N = NGauss;
    M = (N+1)/2;
    x = dvector(1,M);
    w = dvector(1,M);
    gauss_leg_knots();

    /* debug! */
    if (DEBUG) {
      for (i=1; i<=M; i++)
	fprintf(stderr," %d> %25.15f %25.15f\n",i,x[i],w[i]);
    }
  }

/* Do integral */

  bma = 0.5*(b-a);
  bpa = 0.5*(b+a);
  for (i=1, accum=0.0; i<M; i++) {
    accum += w[i]*(*f)( bma*x[i]+bpa);
    accum += w[i]*(*f)(-bma*x[i]+bpa);
  }

  if (N+1-M == M)
    accum += w[M]*(*f)( bma*x[M]+bpa);
  else {
    accum += w[M]*(*f)( bma*x[M]+bpa);
    accum += w[M]*(*f)(-bma*x[M]+bpa);
  }
    

/* Done! */
  return bma*accum;
}

double gint_2(double a, double b, double (*f) (double), int NGauss)
{
  double accum,bma;
  int i;

/* Get values weights and abcissas for Gauss-Legendre integration if needed */

  if (NGauss != N) {
    if (N != 0) {
      free_dvector(x,1,M);
      free_dvector(w,1,M);
    }
    N = NGauss;
    M = (N+1)/2;
    x = dvector(1,M);
    w = dvector(1,M);
    gauss_leg_knots();

    /* debug! */
    if (DEBUG) {
      for (i=1; i<=M; i++)
	fprintf(stderr," %d> %25.15f %25.15f\n",i,x[i],w[i]);
    }
  }

/* Do integral */

  bma = b-a;
  for (i=1, accum=0.0; i<=M; i++)
    accum += 2.0*w[i]*(*f)(a+bma*x[i]*x[i]);

/* Done! */
  return sqrt(b-a)*accum;
}


/* Adapted from Numerical Recipes, Press et al. z*/
    
#define EPS 3.0e-11

void gauss_leg_knots(void)
{
  int j,i;
  double z1,z,pp,p3,p2,p1;
  
  for (i=1;i<=M;i++)  {
    z=cos(3.14159265358979323846*(i-0.25)/(N+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=N;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=N*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    x[i] = z;
    w[i] = 2.0/((1.0-z*z)*pp*pp);
  }
}

#undef EPS
