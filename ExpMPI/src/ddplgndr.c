#include <math.h>

/* Derivative of the associated Legendre polynomial using the recursion
   relation given in G&R */

double plgndr(int, int, double);

double dplgndr(int l, int m, double x)
{
  double tmp,fac;

  if (!l) return 0.0;
  fac = x*x - 1.0;
  if (fac >= 0.0) return 0.0;
  
  tmp = ( plgndr(l+1,m,x)*(l-m+1) - plgndr(l,m,x)*x*(l+1) )/fac;
  
  return tmp;
}
