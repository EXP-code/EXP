#include <math.h>
#include <numerical.h>

#define JMAX 40

double rtbis(func_1d func, double x1, double x2, double xacc)
{
  int j;
  double dx,f,fmid,xmid,rtb;

  f=(*func)(x1);
  fmid=(*func)(x2);
  if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in RTBIS");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5));
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror("Too many bisections in RTBIS");
	
  return 0.0;
}

#undef JMAX

