#include <math.h>
#include <Vector.h>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void get_ultra(int nmax, double l, double x, Vector& p)
{
  double a,b,u1,u2,u;
  int j;

  p[0] = u2 = 1.0;
  p[1] = u  = 2.0*x*(l+1.0);

  for (j=2; j<=nmax; j++) {
    u1 = u2;
    u2 = u;    
    a = 2.0*x*(l+(double)j)/(double)(j);
    b = -(double)(2.0*l + (double)j)/(double)(j);

    p[j] = u = a*u2 + b*u1;
  }

}

