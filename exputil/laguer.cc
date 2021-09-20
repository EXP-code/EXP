#include <complex>
#include <cmath>
#include <Eigen/Eigen>


// const double EPSS=6.0e-8;
const double EPSS=6.0e-16;
const int MAXIT=2000;

void laguer_root(Eigen::VectorXcd& a, int m,
		 std::complex<double> *x, double eps, int polish)
{
  int j,iter;
  double err,dxold,cdx,abx;
  std::complex<double> sq,h,gp,gm,g2,g,b,d,dx,f,x1;

  dxold= fabs(*x);
  for (iter=1; iter<=MAXIT; iter++) {
    b = a[m];
    err = fabs(b);
    d = f = 0.0;
    abx = fabs(*x);
    for (j=m-1; j>=0; j--) {
      f = (*x)*f + d;
      d = (*x)*d + b;
      b = (*x)*b + a[j];
      err = fabs(b) + abx*err;
    }
    err *= EPSS;
    if (fabs(b) <= err) return;
    g = d/b;
    g2 = g*g;
    h = g2 - 2.0*f/b;
    double dm = m;
    sq = sqrt((h*dm - g2)*(dm-1.0));
    gp = g + sq;
    gm = g - sq;
    if (fabs(gp) < fabs(gm)) gp = gm;
    dx = std::complex<double>((double)m,0.0)/gp;
    x1 = *x - dx;
    if (x->real() == x1.real() && x->imag() == x1.imag()) return;
    *x = x1;
    cdx = fabs(dx);
    dxold = cdx;
    if (!polish)
      if (cdx <= eps*fabs(*x)) return;
  }
  puts("Too many iterations in routine LAGUER");
}

