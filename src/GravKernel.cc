#include <cmath>
#include "GravKernel.H"

std::pair<double, double> PlummerSoft::operator()(double r, double eps)
{
  std::pair<double, double> ret;
  ret.first  = std::pow(r*r/(r*r + eps*eps), 1.5);
  ret.second = - std::pow(eps*eps/(r*r + eps*eps), 1.5)/eps;
  if (r > tol*eps) ret.second += -ret.first/r;

  return ret;
}

std::pair<double, double> SplineSoft::operator()(double r, double eps)
{
  std::pair<double, double> ret;
  double x = r/eps;
  if (x<0.5) {
    ret.first = m1(x);
    ret.second = -(fac1 - p1(x))/eps;
    if (x>tol) ret.second += -ret.first/r;
  } else if (x<1.0) {
    ret.first = fac0 + m2(x);
    ret.second = -ret.first/r - (fac2 - p2(x))/eps;
  } else {
    ret.first  = 1.0;
    ret.second = -1.0/r;
  }

  return ret;
}
