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

#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

#include <gaussQ.H>

                              /* if DEBUG is set, weights and abscissas are
                                 printed on stderr */
#ifndef DEBUG
#   define DEBUG 0
#endif

static int N=0,M;             /* N=0 force comp. of knots on first call */
static LegeQuad *lq = NULL;

double gint_0(double a, double b, double (*f) (double), int NGauss)
{
  double accum,bma,bpa;
  int i;

/* Get values weights and abcissas for Gauss-Legendre integration if needed */

  if (NGauss != N) {

    if (lq) delete lq;

    N = NGauss;
    M = (N+1)/2;
    lq = new LegeQuad(M);

    /* debug! */
    if (DEBUG) {
      for (i=1; i<=M; i++)
	cerr << setw(3) << i << "> " 
	     << setw(18) << lq->knot(i)
	     << setw(18) << lq->weight(i)
	     << "\n";
    }
  }

/* Do integral */

  bma = 0.5*(b-a);
  bpa = 0.5*(b+a);
  for (i=1, accum=0.0; i<M; i++) {
    accum += lq->weight(i)*(*f)( bma*lq->knot(i)+bpa);
    accum += lq->weight(i)*(*f)(-bma*lq->knot(i)+bpa);
  }

  if (N+1-M == M)
    accum += lq->weight(M)*(*f)( bma*lq->knot(M)+bpa);
  else {
    accum += lq->weight(M)*(*f)( bma*lq->knot(M)+bpa);
    accum += lq->weight(M)*(*f)(-bma*lq->knot(M)+bpa);
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
    
    if (lq) delete lq;

    N = NGauss;
    M = (N+1)/2;
    lq = new LegeQuad(M);

    /* debug! */
    if (DEBUG) {
      for (i=1; i<=M; i++)
	cerr << setw(3) << i << "> " 
	     << setw(18) << lq->knot(i)
	     << setw(18) << lq->weight(i)
	     << "\n";
    }
  }

/* Do integral */

  bma = b-a;
  for (i=1, accum=0.0; i<=M; i++)
    accum += 2.0*lq->weight(i)*(*f)(a+bma*lq->knot(i)*lq->knot(i));

/* Done! */
  return sqrt(b-a)*accum;
}

