/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine finds the cubic spline coefficients.  The boundary condi-
 *  tions may be one of the following three:
 *       1) "natural," that is, zero second derivatives
 *       2) first derivatives specified
 *       3) third derivatives computed from supplied data
 *
 *
 *  Call sequence:
 *  -------------
 *  void spline(x,y,n,yp1,ypn,y2);
 *
 *  int n;
 *  double x[n],y[n],yp1,ypn,y2[n];
 *
 *  Parameters:
 *  ----------
 *
 *  n        number of supplied grid points
 *  x        abcissa array
 *  y        ordinate array
 *  yp1      boundary condition at j=1
 *  ypn      boundary condition at j=n
 *  y2       array to contain spline coefficients
 *
 *  Returns:
 *  -------
 *
 *  None, spline coefficients returned by pointer
 *
 *  Notes:
 *  -----
 *
 *  If     yp1,yp2 >  1.0e30  boundary conditions (1) natural splines are used
 *  If     yp1,yp2 < -1.0e30  boundary conditions (3) approx. 3rd derivs used
 *  Otherwise                 boundary conditions (2) explicit 2nd derivs used
 *
 *  By:
 *  --
 *  Adopted from Numerical Recipes, Press et al.
 *  Third deriv. boundary condition ---  MDW 11/13/88
 *
 ***************************************************************************/

#ifndef Spline_H
#define Spline_H

#include <cstdlib>
#include <vector>
#include <cmath>

template <class T>
void Spline(const std::vector<T> &x, const std::vector<T> &y, T yp1, T ypn, 
	    std::vector<T> &y2)
{
  int i,k,i1,i2;
  T d1,d2,p,qn,un;
  T sig;
  
  i1 = 0;
  i2 = x.size() - 1;

  std::vector<T> u(i2);
  y2.erase(y2.begin(), y2.end());
  y2 = std::vector<T>(i2);

  //     Boundary conditions obtained by fixing third derivative as computed
  //     by divided differences
  if (yp1 < -0.99e30) {
    y2[i1+0] = 1.0;
    d2  =  ((y[i1+3]-y[i1+2])/(x[i1+3]-x[i1+2]) - (y[i1+2]-y[i1+1])/(x[i1+2]-x[i1+1]))/(x[i1+3]-x[i1+1]);
    d1  =  ((y[i1+2]-y[i1+1])/(x[i1+2]-x[i1+1]) - (y[i1+1]-y[i1+0])/(x[i1+1]-x[i1+0]))/(x[i1+2]-x[i1+0]);
    u[i1+0] = -6.0*(d2-d1)*(x[i1+1]-x[i1+0])/(x[i1+3]-x[i1+0]);
  }
  //     "Normal" zero second derivative boundary conditions */
  else if (yp1 > 0.99e30)
    y2[i1+0] = u[i1+0] = 0.0;
  //      Known first derivative */
  else {
    y2[i1+0] = -0.5;
    u[i1+0] = (3.0/(x[i1+1]-x[i1+0]))*((y[i1+1]-y[i1+0])/(x[i1+1]-x[i1+0])-yp1);
  }
  for (i=i1+1;i<i2;i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  //     Boundary conditions obtained by fixing third derivative as computed
  //     by divided differences
  if (ypn < -0.99e30) {
    d2 = ((y[i2]-y[i2-1])/(x[i2]-x[i2-1]) - 
	  (y[i2-1]-y[i2-2])/(x[i2-1]-x[i2-2]))/(x[i2]-x[i2-2]);
    d1 = ((y[i2-1]-y[i2-2])/(x[i2-1]-x[i2-2]) - 
	  (y[i2-2]-y[i2-3])/(x[i2-2]-x[i2-3]))/(x[i2-1]-x[i2-3]);
    qn = -1.0;
    un = 6.0*(d2-d1)*(x[i2]-x[i2-1])/(x[i2]-x[i2-3]);
  }
  //     "Normal" zero second derivative boundary conditions */
  else if (ypn > 0.99e30)
    qn = un = 0.0;
  //      Known first derivative */
  else {
    qn = 0.5;
    un = (3.0/(x[i2]-x[i2-1]))*(ypn-(y[i2]-y[i2-1])/(x[i2]-x[i2-1]));
  }
  y2[i2] = (un-qn*u[i2-1])/(qn*y2[i2-1]+1.0);
  for (k=i2-1;k>=i1;k--)
    y2[k] = y2[k]*y2[k+1]+u[k];
}



template <class T>
void Splint1(const std::vector<T> &xa, const std::vector<T> &ya, const std::vector<T> &y2a, 
	     T x, T &y, int even=0)
{
  int klo, khi, n1, n2, k;
  T h,b,a;
  
  n1 = 0;
  n2 = xa.size() - 1;

  if (even) {
    klo = (int)( (x-xa[n1])/(xa[n2]-xa[n1])*(double)(n2-n1) ) + n1;
    klo = klo<n1 ? n1 : klo;
    klo = klo<n2 ? klo : n2-1;
    khi = klo+1;
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

  h = xa[khi]-xa[klo];
  
  if (h == 0.0) {
    cerr << "Bad XA input to routine Splint1\n";
    exit(-1);
  }
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

template <class T>
void Splint2(const std::vector<T> &xa, const std::vector<T> &ya, const std::vector<T> &y2a, 
	     T x, T &y, T &dy, int even=0)
{
  int klo, khi, n1, n2, k;
  T h,b,a;

  n1 = 0;
  n2 = xa.size() - 1;

  if (even) {
    klo = (int)( (x-xa[n1])/(xa[n2]-xa[n1])*(double)(n2-n1) ) + n1;
    klo = klo<n2 ? klo : n2-1;
    khi = klo+1;
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

  h = xa[khi]-xa[klo];
  
  if (h == 0.0) {
    cerr << "Bad XA input to routine Splint2\n";
    exit(-1);
  }
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  dy = (-ya[klo]+ya[khi])/h +
    (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
      *h/6.0;
/*  ddy = a*y2a[klo]+b*y2a[khi]; */
}


template <class T>
void Splint3(const std::vector<T> &xa, const std::vector<T> &ya, const std::vector<T> &y2a, 
	     T x, T &y, T &dy, T &ddy, int even=0)
{
  int klo, khi, n1, n2, k;
  T h,b,a;
  
  n1 = 0;
  n2 = xa.size() - 1;

  if (even) {
    klo = (int)( (x-xa[n1])/(xa[n2]-xa[n1])*(double)(n2-n1) ) + n1;
    klo = klo<n2 ? klo : n2-1;
    khi = klo+1;
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

#endif
