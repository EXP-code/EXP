/**************************************************************************
*
*		Computes m zeros of spherical Bessel functions (first kind) 
*               j_n(x) by brute force and using Watson's relation for the
*               location of the first zero
*
*
*		MDW: 12/26/1987
**************************************************************************/
#include <math.h>
#include "cutil.h"

/* test routine 
main()
{
  int n,m,i;
  double *sbessjz(),*a;

  printf("This routine finds m zeros of spherical bessel functions of\n");
  printf("order n . . . n,m: ");
  scanf("%d %d",&n, &m);

  a = sbessjz(n,m);

  for (i=1; i<=m; i++)
    printf(" %d  %.10lf\n",i,a[i]);
}
   end test routine */


#define PI 3.1415926535897932385
#define STEPS 6
#define TOL 1.0e-7

static int NN;

double *sbessjz(int n, int m)
{
  double sbessj(int n, double x);
  double zbrent(double (*func) (double), double x1, double x2, double tol);
  double zbess(double z);
  double *a,z,dz,zl,f,fl;
  int i;

  a = dvector(1,m);

  NN = n;
  dz = PI/STEPS;
  for (i=1, zl=z=0.5+fabs((double)n), fl=sbessj(n,z); i<=m; i++) {
    z += dz;
    f = sbessj(n,z);
    while (f*fl>0) {
      zl = z;
      fl = f;
      z += dz;
      f = sbessj(n,z);
    }
    a[i] = zbrent(zbess,zl,z,TOL);
    zl = z;
    fl = f;
  }

  return a;

}

double zbess(double z)
{
  double sbessj(int n, double x);
  return sbessj(NN,z);
}
