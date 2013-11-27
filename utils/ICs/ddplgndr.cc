#include <math.h>

void nrerror(const char *);


/* Derivative of the associated Legendre polynomial using the recursion
   relation given in G&R */

// This bit from NR
double plgndr2(int l, int m, double x)
{
  double fact, pll, pmm, pmmp1, somx2;
  int i, ll;
  
  if (m < 0 || m > l || fabs(x) > 1.0) {
    nrerror("Bad arguments in routine DPLGNDR2");
  }

  pmm = 1.0;

  if (m > 0) {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact = 1.0;
    for (i=1;i<=m;i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1 = x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=(m+2);ll<=l;ll++) {
	pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      return pll;
    }
  }
}


double dplgndr2(int l, int m, double x)
{
  double tmp,fac;

  if (!l) return 0.0;
  fac = x*x - 1.0;
  if (fac >= 0.0) return 0.0;
  
  tmp = ( plgndr2(l+1,m,x)*(l-m+1) - plgndr2(l,m,x)*x*(l+1) )/fac;
  
  return tmp;
}
