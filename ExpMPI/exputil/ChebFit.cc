//
// Chebyshev fitting and smoothing class
//

#include <iostream>
#include <string>
#include <vector>

#include <math.h>
#include <stdarg.h>
#include <interp.h>
#include <ChebFit.h>

ChebFit::ChebFit(Vector& X, Vector &Y, int N)
{

				// set private values
  n = N;
  a = X[X.getlow()];
  b = X[X.gethigh()];
  rcheck = true;

  new_data(X, Y, N);

}

ChebFit::ChebFit(double A, double B, Vector& X, Vector &Y, int N)
{

				// set private values
  n = N;
  a = A;
  b = B;
  rcheck = true;

  new_data(X, Y, N);

}


void ChebFit::new_limits(double A, double B)
{
  a = A;
  b = B;
}
    

void ChebFit::new_func(double (*func)(double), double A, double B, int N)
{
  double y, sum, fac, bpa, bma;

  a = A;
  b = B;
  n = N;

  Vector f(0,n-1);
  c.setsize(0,n-1);
  c1.setsize(0,n-1);
  c2.setsize(0,n-1);

  bma = 0.5*(b - a);
  bpa = 0.5*(b + a);
  for (int k=0; k<n; k++) {
    y = bma*cos(M_PI*(k+0.5)/n) + bpa;
    f[k] = (*func)(y);
  }

  fac=2.0/n;
  for (int j=0; j<n; j++) {
    sum=0.0;
    for (int k=0;k<n;k++)
      sum += f[k]*cos(M_PI*j*(k+0.5)/n);
    c[j]=fac*sum;
  }

  chder(c, c1);
  chder(c1, c2);

  defined = true;
}


void ChebFit::new_data(Vector& X, Vector& Y, int N)
{
  double y, sum, fac, bpa, bma;

  a = X[X.getlow()];
  b = X[X.gethigh()];
  n = N;

  Vector f(0,n-1);
  c.setsize(0,n-1);
  c1.setsize(0,n-1);
  c2.setsize(0,n-1);

  bma = 0.5*(b - a);
  bpa = 0.5*(b + a);
  for (int k=0; k<n; k++) {
    y = bma*cos(M_PI*(k+0.5)/n) + bpa;
    f[k] = odd2(y, X, Y);
  }

  fac=2.0/n;
  for (int j=0; j<n; j++) {
    sum=0.0;
    for (int k=0; k<n; k++)
      sum += f[k]*cos(M_PI*j*(k+0.5)/n);
    c[j]=fac*sum;
  }

  chder(c, c1);
  chder(c1, c2);

  defined = true;
}


void ChebFit::chder(Vector& cin, Vector& cder)
{
  double con;

  cder[n-1] = 0.0;
  cder[n-2] = 2*(n-1)*cin[n-1];
  for (int j=n-3; j>=0; j--)
    cder[j] = cder[j+2]+ 2*(j+1)*cin[j+1];

  con = 2.0/(b-a);
  for (int j=0; j<n; j++)
    cder[j] *= con;
}


double ChebFit::chebev(double x, Vector& cin)
{
  double d=0.0, dd=0.0, sv, y, y2;
  if (rcheck && (x-a)*(x-b) > 0.0) 
    bomb("ChebFit::chebev", "x out of range", 0);

  y2 = 2.0*(y=(2.0*x-a-b)/(b-a));

  for (int j=n-1;j>=1;j--) {
    sv = d;
    d = y2*d-dd+cin[j];
    dd = sv;
  }

  return y*d-dd+0.5*cin[0];
}

void ChebFit::bomb(const char *a, ...)
{
  va_list ap;
  char *b, *c;

  va_start(ap, a);
  b = va_arg(ap, char *);
  c = va_arg(ap, char *);
  va_end(ap);

  if (b != (char *)0)
    cerr << a << ": " << b;
  if (c != (char *)0)
    cerr << ", " << c;
  cerr << '\n';

  exit(-1);
}
