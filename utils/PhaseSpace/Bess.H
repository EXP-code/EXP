#ifndef _Bess_H
#define _Bess_H

#include <cmath>

#include <EXPmath.H>

//! Generate orthonormal Bessel functions of integral order
class Bess
{
private:

  unsigned int nroots;
  double order;
  std::vector<double> roots, norm;

  std::tuple<double, double>
  bracket(std::function<double(double)> B, double zmin, double z0);

  std::vector<double>
  zeros(double Nu, int Nz, double tol=1.0e-18);

public:

  //! Constructor: set the order and number of radial functions
  Bess(double order, unsigned int nroots) : order(order), nroots(nroots)
  {
    roots = zeros(order, nroots);
    norm.resize(nroots);
    for (unsigned int m=0; m<nroots; m++) {
      double val = EXPmath::cyl_bessel_j(order+1.0, roots[m]);
      norm[m] = sqrt(0.5*val*val);
    }
  }
  
  //! Get the norm for radial order m
  double getNorm(int m)
  {
    if (m>=nroots) return 0.0;
    else           return norm[m];
  }
  
  //! Evaluate the Bessel for x in [0, 1] for radial order m
  double operator()(double& x, const unsigned& m)
  {
    if (m>=nroots) return 0.0;
    return EXPmath::cyl_bessel_j(order, x*roots[m]) / norm[m];
  } 

  //! Evaluate the Bessel for x in [0, 1] for radial order m
  double eval(double& x, unsigned& m)
  {
    if (m>=nroots) return 0.0;
    return EXPmath::cyl_bessel_j(order, x*roots[m]) / norm[m];
  } 

}; 

#endif
