/**************************************************************************
*
*		Computes m zeros of spherical Bessel functions (first kind) 
*               j_n(x) by brute force and using Watson's relation for the
*               location of the first zero
*
*
*		MDW: 12/26/1987
*                    port to C++ 02/15/94
**************************************************************************/

#include <cmath>
#include <numerical.H>
#include <EXPmath.H>

#define STEPS 6
#define TOL 1.0e-7

static int NN;

static double zbess(double z)
{
  return EXPmath::sph_bessel(NN, z);
}

Eigen::VectorXd sbessjz(int n, int m)
{
  Eigen::VectorXd a(m);

  auto zfunc = [n](double z) { return std::cyl_bessel_j(n, z); };

  double dz = M_PI/STEPS;
  double z  = 0.5+fabs((double)n);
  double zl = z, fl, f;
  for (int i=0; i<m; i++) {
    fl = EXPmath::sph_bessel(n,z);
    z += dz;
    f  = EXPmath::sph_bessel(n, z);
    while (f*fl>0) {
      zl = z;
      fl = f;
      z += dz;
      f = EXPmath::sph_bessel(n,z);
    }
    a[i] = zbrent(zbess, zl, z, TOL);
    zl = z;
    fl = f;
  }

  return a;
}
