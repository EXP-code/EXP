
#include <math.h>
#include <iostream.h>
#include <stdlib.h>
#include <Vector.h>

void Splint1(const Vector &xa, const Vector &ya, const Vector &y2a, double x, double &y, int even=0)
{
  int klo, khi, n1, n2, k;
  double h,b,a;
  void nrerror(char *);
  
  n1 = xa.getlow();
  n2 = xa.gethigh();

  if (even) {
    klo=(int)( (x-xa[n1])/(xa[n2]-xa[n1])*(double)(n2-n1) ) + n1;
    klo=klo<n1 ? n1 : klo;
    klo=klo<n2 ? klo : n2-1;
    khi=klo+1;
  }
  else {
    klo = n1;
    khi = n2;
    while (khi-klo > 1) {
      k = (khi+klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
  }

  h=xa[khi]-xa[klo];
  
  if (h == 0.0) {
    cerr << "Bad XA input to routine Splint1\n";
    exit(-1);
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
/*  yd=(-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
      *h/6.0;
  ydd = a*y2a[klo]+b*y2a[khi]; */
}

void Splint2(const Vector &xa, const Vector &ya, const Vector &y2a, double x, double &y, double &dy, int even=0)
{
  int klo, khi, n1, n2, k;
  double h,b,a;
  void nrerror(char *);
  
  n1 = xa.getlow();
  n2 = xa.gethigh();

  if (even) {
    klo=(int)( (x-xa[n1])/(xa[n2]-xa[n1])*(double)(n2-n1) ) + n1;
    klo=klo<n2 ? klo : n2-1;
    khi=klo+1;
  }
  else {
    klo = n1;
    khi = n2;
    while (khi-klo > 1) {
      k = (khi+klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
  }

  h=xa[khi]-xa[klo];
  
  if (h == 0.0) {
    cerr << "Bad XA input to routine Splint2\n";
    exit(-1);
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  dy=(-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
      *h/6.0;
/*  ddy = a*y2a[klo]+b*y2a[khi]; */
}


void Splint3(const Vector &xa, const Vector &ya, const Vector &y2a, double x, double &y, double &dy, double &ddy, int even=0)
{
  int klo, khi, n1, n2, k;
  double h,b,a;
  void nrerror(char *);
  
  n1 = xa.getlow();
  n2 = xa.gethigh();

  if (even) {
    klo=(int)( (x-xa[n1])/(xa[n2]-xa[n1])*(double)(n2-n1) ) + n1;
    klo=klo<n2 ? klo : n2-1;
    khi=klo+1;
  }
  else {
    klo = n1;
    khi = n2;
    while (khi-klo > 1) {
      k = (khi+klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
  }

  h=xa[khi]-xa[klo];
  
  if (h == 0.0) {
    cerr << "Bad XA input to routine Splint3\n";
    exit(-1);
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  dy=(-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
      *h/6.0;
  ddy = a*y2a[klo]+b*y2a[khi];
}
