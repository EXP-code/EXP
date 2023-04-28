#ifndef _HankelTransform_H
#define _HankelTransform_H

#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/tools/minima.hpp>

//! Hankel transform based on the Ogata quadrature formulae from:
//!
//! Ogata, H, 2005, "A Numerical Integration Formula Based on the
//! Bessel Functions", Publ. Res. Inst. Math. Sci. 41 (4) (2005)
//! 949–970.
//!
class HankelTransform
{
private:
  //! the h scaling value for the Ogata mapping
  double h;

  //! nu is Bessel function order
  double nu;

  //! N is number of function calls
  int N;

  //! Will contain the zeros of the Bessel function of order nu
  std::vector<double> zeros;

  //! Zeros divided by 2*pi for Ogata
  std::vector<double> xi;

  //! Ogata weights; Bessel functions of the second kind divided by
  //! the first kind for nu+1 evaluated at the zeros
  std::vector<double> w;

  //! True for mapping transformation algorithm (default)
  bool mapped;

  //! The mapping transformation algorithm (default)
  double ogata_transformed
  (std::function<double (double) > f, double q, double h);

  //! The unmapped algorithm
  double ogata_linear
  (std::function<double (double) > f, double q, double h);

  //! Variable change
  double fk_trans(double x, std::function<double (double) > f, double q)
  {
    return f(x/q)/q;
  };

  //! Ogata mapping formula
  double get_psi(double t);

  //! Derivative of the Ogata mapping formula
  double get_psip(double t);

public:
  //! Null constructor
  HankelTransform() {}

  //! Constructor
  HankelTransform(double h, double nu, int N=20);

  //! Destructor
  ~HankelTransform() {}

  //! Perform the transform F(q)=int(f(x)*Jn(x*q))
  double operator()(std::function<double (double) > f,  double q);

  //! Set to unmapped Ogata formula
  void setLinear() { mapped = false; }

};

#endif
// END: _HankelTransform_H
