#ifndef _HankelTransform_H
#define _HankelTransform_H

#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/tools/minima.hpp>

class HankelTransform
{
private:
  // the h scaling value for the Ogata mapping
  double h;

  // nu is Bessel function order
  double nu;

  // N is number of function calls
  int N;

  std::vector<double> jn_zeros0;
  std::vector<double> zeros;
  std::vector<double> xi;
  std::vector<double> Jp1;
  std::vector<double> w;

  bool mapped;

  double ogata_transformed
  (std::function<double (double) > f, double q, double h);

  double ogata_linear
  (std::function<double (double) > f, double q, double h);

  double fk_trans(double x, std::function<double (double) > f, double q)
  {
    return f(x/q)/q;
  };

  double get_psi(double t);
  double get_psip(double t);

  bool transformed;

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
