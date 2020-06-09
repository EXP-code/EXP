#ifndef _GravKernel_H
#define _GravKernel_H

#include <utility>

class PlummerSoft
{
public:
  std::pair<double, double> operator()(double r, double eps);

};

class SplineSoft
{
private:
  double eps, fac1, fac2;

  double m1(double x)
  { return 32.*x*x*x*(1./3. - 6./5.*x*x + x*x*x); }

  double m2(double x)
  { return 16./15.*x*x*x*(20. - 45.*x + 36.*x*x - 10.*x*x*x); }

  double p1(double x)
  { return 32.*x*x*(0.5 - 1.5*x*x + 6./5.*x*x*x); }

  double p2(double x)
  { return 32.*x*x*(1. - 2.*x + 1.5*x*x - 2./5.*x*x*x); }

public:
  SplineSoft()
  {
    fac1 = p2(1.0) - p2(0.5) + p1(0.5);
    fac2 = p2(1.0);
  }

  std::pair<double, double> operator()(double r, double eps);
};

#endif
