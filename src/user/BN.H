#ifndef _BN_H
#define _BN_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

/**
   Cumulative angle scattering from Bobylev and Nanbu 2000
*/
class coulombSelect
{
private:

  double f(double x)
  {
    return 1.0/tanh(x) - 1.0/x;
  }

  double deriv(double x)
  {
    if (x>100.0) return 1.0/(x*x);
    double z = sinh(x);
    return 1.0/(x*x) - 1.0/(z*z);
  }

  double func(double tau, double x)
  {
    const double tol = 1.0e-12;
    const int itmax  = 1000;
    double fac       = exp(-2.0*tau), corr;

    for (int i=0; i<itmax; i++) {
      corr  = (f(x) - fac)/deriv(x);
      x    -= corr;
      if (fabs(corr)<tol) break;
    }

    return x;
  }

  const double tau_i = 0.0001, tau_f = 4.0, tau_z = 40.0;
  const int numT = 2000;
  std::vector<double> cT, cA;

public:
  
  //! Constructor
  coulombSelect()
  {
    cT.resize(numT);
    cA.resize(numT);

    double del = (log(tau_f) - log(tau_i))/(numT-1);
    double A   = 0.5/tau_i;
    for (int i=0; i<numT; i++) {
      cT[i] = tau_i*exp(del*i);
      cA[i] = A = func(cT[i], A);
    }
  }

  //! Evaluate A given tau
  double eval(double tau)
  {
    if      (tau <= tau_i) return 1.0/(2.0*tau);
    else if (tau >= tau_z) return 0.0;
    else if (tau >= tau_f) return 3.0*exp(-2.0*tau);
    else {
      auto it    = std::lower_bound(cT.begin(), cT.end(), tau);
      int indx   = std::distance(cT.begin(), it);

      double lo  = cT[indx-1];
      double hi  = cT[indx];
      double del = hi - lo;
      double A   = (hi - tau)/del;
      double B   = (tau - lo)/del;

      return A*cA[indx-1] + B*cA[indx];
    }
  }

  //! Select tau given random number U in [0,1)
  double operator()(double tau, double U)
  {
    double A = eval(tau);
    if (A<1.0e-10)
      return 2.0*U - 1.0;
    else if (A>40.0)
      return 1.0 + log(U)/A;
    else
      return log(exp(-A) + 2.0*U*sinh(A))/A;
  }

};
  
#endif // _BN_H

