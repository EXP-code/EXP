#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Eigen/Eigen>

using namespace std;

void Splint1(const Eigen::VectorXd &xa,
	     const Eigen::VectorXd &ya,
	     const Eigen::VectorXd &y2a, double x, double &y)
{
  int sz  = xa.size();
  int k   = 0;
  int klo = 0;
  int khi = sz-1;

  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }

  double h = xa[khi] - xa[klo];
  
  if (h == 0.0) {
    std::cerr << "Bad XA input to routine Splint1\n";
    exit(-1);
  }
  double a = (xa[khi]-x)/h;
  double b = (x-xa[klo])/h;

  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void Splint2(const Eigen::VectorXd &xa,
	     const Eigen::VectorXd &ya,
	     const Eigen::VectorXd &y2a, double x, double &y, double &dy)
{
  int sz  = xa.size();
  int k   = 0;
  int klo = 0;
  int khi = sz-1;

  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }

  double h = xa[khi] - xa[klo];
  
  if (h == 0.0) {
    cerr << "Bad XA input to routine Splint2\n";
    exit(-1);
  }

  double a=(xa[khi]-x)/h;
  double b=(x-xa[klo])/h;

  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

  dy = (-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])*h/6.0;
}


void Splint3(const Eigen::VectorXd &xa,
	     const Eigen::VectorXd &ya,
	     const Eigen::VectorXd &y2a,
	     double x, double &y, double &dy, double &ddy)
{
  int sz = xa.size();

  int klo = 0;
  int khi = sz-1;
  int k   = 0;

  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }

  double h = xa[khi] - xa[klo];
  
  if (h == 0.0) {
    std::cerr << "Bad XA input to routine Splint3\n";
    exit(-1);
  }
  double a = (xa[khi]-x)/h;
  double b = (x-xa[klo])/h;

  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

  dy = (-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])*h/6.0;

  ddy = a*y2a[klo]+b*y2a[khi];
}
