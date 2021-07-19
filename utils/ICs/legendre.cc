#include <iostream>
#include <cmath>

#include <Eigen/Eigen>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

void legendre_R(int lmax, double x, Eigen::MatrixXd& p)
{
  double fact, somx2, pll, pl1, pl2;

  p(0, 0) = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (int m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p(m, m) = pll;
      if (std::isnan(p(m, m)))
	std::cerr << "legendre_R: p[" << m << "][" << m << "]: pll="
		  << pll << std::endl;
      fact += 2.0;
    }
  }

  for (int m=0; m<lmax; m++) {
    pl2 = p(m, m);
    p(m+1, m) = pl1 = x*(2*m+1)*pl2;
    for (int l=m+2; l<=lmax; l++) {
      p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (std::isnan(p(l, m)))
	std::cerr << "legendre_R: p[" << l << "][" << m << "]: pll="
		  << pll << std::endl;

      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (std::isnan(x))
    std::cerr << "legendre_R: x\n";

  for(int l=0; l<=lmax; l++)
    for (int m=0; m<=l; m++)
      if (std::isnan(p(l, m)))
	std::cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
		  << lmax << std::endl;

}

void dlegendre_R(int lmax, double x, Eigen::MatrixXd &p, Eigen::MatrixXd &dp)
{
  const double tol = 1.0e-8;
  double fact, somx2, pll, pl1, pl2;

  if (x < -1.0+tol) x = -1.0+tol;
  if (x >  1.0-tol) x =  1.0-tol;

  p(0, 0) = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (int m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p(m, m) = pll;
      fact += 2.0;
    }
  }

  for (int m=0; m<lmax; m++) {
    pl2 = p(m, m);
    p(m+1, m) = pl1 = x*(2*m+1)*pl2;
    for (int l=m+2; l<=lmax; l++) {
      p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }

  somx2 = 1.0/(x*x - 1.0);
  dp(0, 0) = 0.0;
  for (int l=1; l<=lmax; l++) {
    for (int m=0; m<l; m++)
      dp(l, m) = somx2*(x*l*p(l, m) - (l+m)*p(l-1, m));
    dp(l, l) = somx2*x*l*p(l, l);
  }
}
