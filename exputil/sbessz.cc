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

#define STEPS 6
#define TOL 1.0e-7

static int NN;

static double zbess(double z)
{
  return std::sph_bessel(NN, z);
}

Eigen::VectorXd sbessjz(int n, int m)
{
  double z,dz,zl,f,fl;
  int i;

  Eigen::VectorXd a(m);

  NN = n;
  dz = M_PI/STEPS;
  for (int i=0, zl=z=0.5+fabs((double)n), fl=std::sph_bessel(n,z); i<m; i++) {
    z += dz;
    f = std::sph_bessel(n, z);
    while (f*fl>0) {
      zl = z;
      fl = f;
      z += dz;
      f = std::sph_bessel(n,z);
    }
    a[i] = zbrent(zbess, zl, z, TOL);
    zl = z;
    fl = f;
  }

  return a;
}
