#include <stdio.h>
#include <math.h>
#include "cutil.h"

void usage(char *dummy){};

double ultra(int, double, double);
void get_ultra(int, double, double, double *);

main()
{
  double *u, *du, x, tmp, del;
  int nmax, l, i;

  del = 0.0001;
  x = 0.35;
  l = 4;
  nmax = 12;

  u = dvector(0,nmax);
  du = dvector(0,nmax);

  get_ultra(nmax-1, (double)l,   x, u);
  get_ultra(nmax-1, (double)(l+1),   x, du);
  for (i=nmax; i>0; i--) du[i] = du[i-1];
  du[0] = 0.0;

  for (i=1; i<=nmax; i++) {
    tmp = (ultra(i-1, l, x+del) - ultra(i-1, l, x-del))/(2.0*del);
    printf(" %d>  %e  %e\n", i, tmp, 2.0*(l+1)*du[i-1]);
  }

}
