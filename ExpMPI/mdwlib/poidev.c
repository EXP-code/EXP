#define _BSD_SOURCE
#include <math.h>

double poidev(double xm)
{
  static double sq,alxm,g,oldm=(-1.0);
  double em,t,y;
  double rnd_01d(void),dgammln(double xx);
  
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    t=1.0;
    do {
      em += 1.0;
      t *= rnd_01d();
    } while (t > g);
  } else {
    if (xm != oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-dgammln(xm+1.0);
    }
    do {
      do {
	y=tan(M_PI*rnd_01d());
	em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-dgammln(em+1.0)-g);
    } while (rnd_01d() > t);
  }
  return em;
}
