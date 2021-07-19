#include <complex>
#include <cmath>
#include <Eigen/Eigen>

void laguer_root(Eigen::VectorXcd& a, int m,
		 std::complex<double> *x, double eps, int polish);

double EPS=2.0e-16;

void zroots(Eigen::VectorXcd& a, Eigen::VectorXcd& roots, int polish)
{
  int jj, j, i;
  int m = roots.size();
  std::complex<double> x, b, c;
  
  Eigen::VectorXcd ad = a;
  for (int j=m-1; j>=0; j--) {
    x = 0.0;
    laguer_root(ad, j+1, &x, EPS, 0);
    if (fabs(x.imag()) <= (2.0*EPS*fabs(x.real()))) x.imag(0.0);
    roots[j]=x;
    b=ad[j];
    for (int jj=j-1; jj>=0; jj--) {
      c=ad[jj];
      ad[jj]=b;
      b = x*b + c;
    }
  }
  if (polish)
    for (int j=0;j<m;j++)
      laguer_root(a, m, &roots[j], EPS, 1);
  for (int j=1; j<m; j++) {
    x=roots[j];
    for (int i=j-1;i>=0; i--) {
      if (roots[i].real() <= x.real()) break;
      roots[i+1]=roots[i];
    }
    roots[i+1]=x;
  }
}

