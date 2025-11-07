#include "OrthoPoly.H"

double OrthoPoly::f(const double x, const int n)
{
  if (n==0) return f0(x);
  if (n==1) return f1(x);

  double g2;
  double g1 = f0(x);
  double g0 = f1(x);
  
  for (int nn=2; nn<=n; nn++) {
    g2 = g1;
    g1 = g0;

    g0 = ( (coef2(nn-1) + coef3(nn-1)*x)*g1 - coef4(nn-1)*g2 )/
      coef1(nn-1);
  }

  return g0;
}

Eigen::VectorXd OrthoPoly::fv(const double x, const int n)
{
  Eigen::VectorXd t(n+1);

  t[0] = f0(x);
  t[1] = f1(x);

  double g2;
  double g1 = t[0];
  double g0 = t[1];
  
  for (int nn=2; nn<=n; nn++) {
    g2 = g1;
    g1 = g0;

    g0 = ( (coef2(nn-1) + coef3(nn-1)*x)*g1 - coef4(nn-1)*g2 )/
      coef1(nn-1);
    
    t[nn] = g0;
  }

  return t;
}

