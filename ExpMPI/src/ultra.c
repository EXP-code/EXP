#include <math.h>

#ifdef RCSID
static char rcsid[] = "$Id";
#endif

double ultra(int n, double l, double x)
{
  double a,b,u1,u2,u;
  int j;


  switch (n) {
  case 0:
    return 1.0;
  case 1:
    return 2.0*x*(l+1.0);
  }

  u2 = 1.0;
  u  = 2.0*x*(l+1.0);

  for (j=2; j<=n; j++) {
    u1 = u2;
    u2 = u;    
    a = 2.0*x*(l+(double)j)/(double)(j);
    b = -(double)(2.0*l + (double)j)/(double)(j);

    u = a*u2 + b*u1;
  }

  return u;

}

