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

#include <math.h>
#include <Vector.h>
#include <numerical.h>

#define STEPS 40
#define TOL 1.0e-10

static int NN;

static double zbess(double z)
{
  return jn(NN, z)*NN/z - jn(NN+1, z);
}

Vector bessjz(int n, int m)
{
  double z,dz,zl,f,fl;
  int i;

  Vector a(1, m);

  NN = n;
  dz = M_PI/STEPS;
  for (i=1, zl=z=1.0e-8, fl=zbess(z); i<=m; i++) {
    z += dz;
    f = zbess(z);
    while (f*fl>0) {
      zl = z;
      fl = f;
      z += dz;
      f = zbess(z);
    }
    a[i] = zbrent(zbess,zl,z,TOL);
    z = a[i];
    zl = a[i] + 100.0*TOL;
    fl = zbess(zl);
  }

  return a;

}

