#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <math.h>
  
void sinecosine_R(int mmax, double phi, double *c, double *s)
{
  int m;

  c[0] = 1.0;
  s[0] = 0.0;

  c[1] = cos(phi);
  s[1] = sin(phi);

  for (m=2; m<=mmax; m++) {
    c[m] = 2.0*c[1]*c[m-1] - c[m-2];
    s[m] = 2.0*c[1]*s[m-1] - s[m-2];
  }
}

void sinecosine(int mmax, double phi, double **cc, double **ss)
{
  double *s, *c;
  int m;

  c = *cc;
  s = *ss;

  c[0] = 1.0;
  s[0] = 0.0;

  c[1] = cos(phi);
  s[1] = sin(phi);

  for (m=2; m<=mmax; m++) {
    c[m] = 2.0*c[1]*c[m-1] - c[m-2];
    s[m] = 2.0*c[1]*s[m-1] - s[m-2];
  }
}


