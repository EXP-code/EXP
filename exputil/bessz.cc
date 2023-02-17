/**************************************************************************
*
*		Computes m zeros of spherical Bessel functions (first
*               kind) j_n(x) by brute force and using Watson's
*               relation for the location of the first zero
*
*
*		MDW: 12/26/1987
*                    port to C++ 02/15/94
*                    updated for C++17 01/11/22
**************************************************************************/

#include <cmath>
#include <numerical.H>
#include <EXPmath.H>

#define STEPS 6
#define TOL 1.0e-7

Eigen::VectorXd bessjz(int n, int m)
{
  Eigen::VectorXd a(m);

  auto zfunc = [n](double z) { return EXPmath::cyl_bessel_j(n, z); };

  double dz = M_PI/STEPS;
  double z  = 0.5+fabs((double)n), zl, fl, f;

  for (int i=0; i<m; i++) {
    fl = EXPmath::cyl_bessel_j(n, z);
    z += dz;
    f  = EXPmath::cyl_bessel_j(n, z);
    while (f*fl>0.0) {
      zl = z;
      fl = f;
      z += dz;
      f = EXPmath::cyl_bessel_j(n,z);
    }
    a[i] = zbrent(zfunc, zl, z, TOL);
    zl = z;
    fl = f;
  }

  return a;
}
