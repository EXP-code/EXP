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

#define STEPS 6
#define TOL 1.0e-7

static int NN;

double sbessj(int n, double x);
double zbrent(double (*func) (double), double x1, double x2, double tol);

double zbess(double z)
{
  return sbessj(NN,z);
}

double *sbessjz(int n, int m)
{
  double *a,z,dz,zl,f,fl;
  int i;

  a = new double [m];
  a -= 1;

  NN = n;
  dz = M_PI/STEPS;
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

/* test routine 
int main()
{
  int n,m,i;
  double *a;

  cout << "This routine finds m zeros of spherical bessel functions of\n";
  cout << "order n . . . n,m: ";
  cin >> n;
  cin >> m;

  a = sbessjz(n,m);

  for (i=1; i<=m; i++)
  cout << setw(5) << m << setw(15) << a[i] << "\n";

  delete (a+1);
}
end test routine */

