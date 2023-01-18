#include <math.h>

using namespace std;

#include <functional>
#include <iostream>
#include <iomanip>
#include <vector>

enum ZBrentReturn {
  Good, 
  Bracket, 
  Iteration
};

template <class T>
class ZBrent 
{
  using func_T = std::function<double(double, T&)>;

public:
  static int ITMAX;
  static double EPS;
  static double epsilon;

  ZBrent() {}
  ZBrentReturn find(func_T func, T& param, double x1, double x2, double tol,
		    double& result);
};


template<class T> int ZBrent<T>::ITMAX=100;
template<class T> double ZBrent<T>::EPS=1.0e-16;
template<class T> double ZBrent<T>::epsilon=1.0e-3;

template <class T>
ZBrentReturn
ZBrent<T>::find(func_T func, T& param, double x1, double x2, double tol,
		double& result)
{
  double a=x1, b=x2, c=0.0, d=0.0, e=0.0, min1, min2;
  double fc, p, q, r, s, tol1, xm;
  double fa=func(a, param);
  double fb=func(b, param);

  if (fb*fa > 0.0) {
    if (a > b) {
      a *= (1+epsilon);
      b *= (1-epsilon);
    } else {
      a *= (1-epsilon);
      b *= (1+epsilon);
    }      
    fa=func(a, param);
    fb=func(b, param);
    if (fb*fa > 0.0) return Bracket;
  }

  fc = fb;
  c = b;

  for (int iter=1; iter<=ITMAX; iter++) {
    if (fb*fc > 0.0) {
      c  = a;
      fc = fa;
      e  = d = b-a;
     }
    if (fabs(fc) < fabs(fb)) {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0*EPS*fabs(b)+0.5*tol;
    xm = 0.5*(c - b);

    if (fabs(xm) <= tol1 || fb == 0.0) {
      result =  b;
      return Good;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb/fa;
      if (a == c) {
	p = 2.0*xm*s;
	q = 1.0 - s;
      } else {
	q = fa/fc;
	r = fb/fc;
	p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p = fabs(p);
      min1 = 3.0*xm*q-fabs(tol1*q);
      min2 = fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e = d;
	d = p/q;
      } else {
	d = xm;
	e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;

    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));

    fb = func(b, param);
  }

  return Iteration;

}

