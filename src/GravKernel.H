#ifndef _GravKernel_H
#define _GravKernel_H

#include <utility>

//! Abstract class for smoothing kernel.  Derive all new kernels from this class.
class SoftKernel
{
public:
  //! Separations smaller than this value are considered to be the same particle
  const double tol = 1.0e-8;

  //! This operator returns the fractional mass and gravitational
  //! potential inside of radius @param r for softening @param eps
  virtual std::pair<double, double> operator()(double r, double eps) = 0;

};


//! Plummer softened gravity (infinite support)
class PlummerSoft :  public SoftKernel
{
public:

  //! Main operator returning enclosed mass and gravitational
  //! potential
  std::pair<double, double> operator()(double r, double eps);

};

//! Cubic-spline softened gravity (compact support)
class SplineSoft : public SoftKernel
{
private:
  double eps, fac0, fac1, fac2;

  //@{
  //! Spline kernel integrals
  
  double m1(double x)
  { return 32.*x*x*x*(1./3. - 6./5.*x*x + x*x*x); }

  double m2(double x)
  { return 16./15.*x*x*x*(20. - 45.*x + 36.*x*x - 10.*x*x*x); }

  double p1(double x)
  { return 32.*x*x*(0.5 - 1.5*x*x + 6./5.*x*x*x); }

  double p2(double x)
  { return 32.*x*x*(1. - 2.*x + 1.5*x*x - 2./5.*x*x*x); }
  //@}

public:

  //! Constructor
  SplineSoft()
  {
    fac0 = m1(0.5) - m2(0.5);
    fac1 = p2(1.0) - p2(0.5) + p1(0.5);
    fac2 = p2(1.0);
  }

  //! Main operator returning enclosed mass and gravitational
  //! potential
  std::pair<double, double> operator()(double r, double eps);
};

#endif
