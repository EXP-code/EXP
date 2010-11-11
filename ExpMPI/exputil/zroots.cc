#include <stdio.h>
#include <math.h>
#include "kevin_complex.h"

void laguer_root(CVector& a, int m, KComplex *x, double eps, int polish);

double EPS=2.0e-16;

void zroots(CVector& a, CVector& roots, int polish)
{
  int jj, j, i;
  int m = roots.gethigh();
  KComplex x, b, c;
  
  CVector ad = a;
  for (j=m; j>=1; j--) {
    x = 0.0;
    laguer_root(ad, j, &x, EPS, 0);
    if (fabs(x.imag()) <= (2.0*EPS*fabs(x.real()))) x.imag()=0.0;
    roots[j]=x;
    b=ad[j];
    for (jj=j-1;jj>=0;jj--) {
      c=ad[jj];
      ad[jj]=b;
      b = x*b + c;
    }
  }
  if (polish)
    for (j=1;j<=m;j++)
      laguer_root(a, m, &roots[j], EPS, 1);
  for (j=2;j<=m;j++) {
    x=roots[j];
    for (i=j-1;i>=1;i--) {
      if (roots[i].real() <= x.real()) break;
      roots[i+1]=roots[i];
    }
    roots[i+1]=x;
  }
}

