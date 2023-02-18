#include <cmath>
#include <iostream>
#include <cstdlib>
#include <Eigen/Eigen>

using namespace std;

void Splint1(const Eigen::VectorXd &xa,
	     const Eigen::VectorXd &ya,
	     const Eigen::VectorXd &y2a, double x, double &y, int even=0)
{
  int klo, khi, n1, n2, k;
  double h,b,a;
  
  int sz = xa.size();

  if (even) {
    klo=(int)( (x-xa[0])/(xa[sz-1]-xa[0])*(double)(sz-1) );
    klo=klo<0 ? 0 : klo;
    klo=klo<sz-1 ? klo : sz-2;
    khi=klo+1;
  }
  else {
    klo = 0;
    khi = sz-1;
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

  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void Splint2(const Eigen::VectorXd &xa,
	     const Eigen::VectorXd &ya,
	     const Eigen::VectorXd &y2a, double x, double &y, double &dy, int even=0)
{
  int klo, khi, k;
  double h, b, a;
  
  int sz = xa.size();

  if (even) {
    klo=(int)( (x-xa[0])/(xa[sz-1]-xa[0])*(double)(sz-1) );
    klo=klo<0 ? 0 : klo;
    klo=klo<sz-1 ? klo : sz-2;
    khi=klo+1;
  }
  else {
    klo = 0;
    khi = sz-1;
    while (khi-klo > 1) {
      k = (khi+klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
  }

  h=xa[khi]-xa[klo];
  
  if (h == 0.0) {
    cerr << "klo=" << klo << " khi=" << khi << " (lo, hi)=(" << xa[klo] << ", " << xa[khi] << ") sz=" << sz << " xa[0]=" << xa[0] << " xa[1]" << xa[1] << std::endl;
    cerr << "Bad XA input to routine Splint2\n";
    exit(-1);
  }
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  dy = (-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
    *h/6.0;
}


void Splint3(const Eigen::VectorXd &xa,
	     const Eigen::VectorXd &ya,
	     const Eigen::VectorXd &y2a,
	     double x, double &y, double &dy, double &ddy, int even=0)
{
  int klo, khi, k;
  double h,b,a;
  
  int sz = xa.size();

  if (even) {
    klo = (int)( (x-xa[0])/(xa[sz-1]-xa[0])*(double)(sz-1) );
    klo = klo<0 ? 0 : klo;
    klo = klo<sz-1 ? klo : sz-2;
    khi = klo+1;
  }
  else {
    klo = 0;
    khi = sz-1;
    while (khi-klo > 1) {
      k = (khi+klo) >> 1;
      if (xa[k] > x) khi = k;
      else klo = k;
    }
  }

  h = xa[khi]-xa[klo];
  
  if (h == 0.0) {
    cerr << "Bad XA input to routine Splint3\n";
    exit(-1);
  }
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  dy = (-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
      *h/6.0;
  ddy = a*y2a[klo]+b*y2a[khi];
}
