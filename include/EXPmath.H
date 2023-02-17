// Work around for incomplete C++17 special function implementations.
// Yes, CLANG, I mean you!

#ifndef _EXPmath_H
#define _EXPmath_H

namespace EXPmath
{
  double cyl_bessel_i(double nu,  double x);
  double cyl_bessel_j(double nu,  double x);
  double cyl_bessel_k(double nu,  double x);
  double sph_bessel  (unsigned n, double x);
}
#endif
