#include <math.h>

#define ITMAX 100
#define EPS 3.0e-7

void gcf(double *gammcf, double a, double x, double *gln)
{
  int n;
  double gold=0.0,g,fac=1.0,b1=1.0;
  double b0=0.0,anf,ana,an,a1,a0=1.0;
  double dgammln(double);
  void nrerror();

  *gln=dgammln(a);
  a1=x;
  for (n=1;n<=ITMAX;n++) {
    an=(double) n;
    ana=an-a;
    a0=(a1+a0*ana)*fac;
    b0=(b1+b0*ana)*fac;
    anf=an*fac;
    a1=x*a0+anf*a1;
    b1=x*b0+anf*b1;
    if (a1) {
      fac=1.0/a1;
      g=b1*fac;
      if (fabs((g-gold)/g) < EPS) {
	*gammcf=exp(-x+a*log(x)-(*gln))*g;
	return;
      }
      gold=g;
    }
  }
  nrerror("a too large, ITMAX too small in routine GCF");
}

#undef ITMAX
#undef EPS
