// Work around for incomplete C++17 special function implementations.

// Yes, CLANG, I mean you!  Come on...

// If you are using anything but CLANG, you will get the C++17 STL
// implementation of special functions.

#ifndef _EXPmath_H
#define _EXPmath_H

namespace AltMath
{
  double cyl_bessel_i  (double nu,  double x);
  double cyl_bessel_j  (double nu,  double x);
  double cyl_bessel_k  (double nu,  double x);
  double sph_bessel    (unsigned n, double x);
  double assoc_laguerre(unsigned l, unsigned n, double x);
}

#ifdef __clang__
namespace EXPmath = AltMath;
#else
namespace EXPmath = std;
#endif


#endif
