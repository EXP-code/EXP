#include <iostream>
#include <math.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <Vector.h>
  
void legendre_R(int lmax, double x, Matrix& p)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p[0][0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[m][m] = pll;
      if (isnan(p[m][m]))
	cerr << "legendre_R: p[" << m << "][" << m << "]: pll=" << pll << "\n";
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "]: pll=" << pll << "\n";

      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (isnan(x))
    cerr << "legendre_R: x\n";
  for(l=0; l<=lmax; l++)
    for (m=0; m<=l; m++)
      if (isnan(p[l][m]))
	cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
	     << lmax << "\n";

}

void dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp)
{
  const double tol = 1.0e-8;
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  if (x < -1.0+tol) x = -1.0+tol;
  if (x >  1.0-tol) x =  1.0-tol;

  p[0][0] = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p[m][m] = pll;
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }

  somx2 = 1.0/(x*x - 1.0);
  dp[0][0] = 0.0;
  for (l=1; l<=lmax; l++) {
    for (m=0; m<l; m++)
      dp[l][m] = somx2*(x*l*p[l][m] - (l+m)*p[l-1][m]);
    dp[l][l] = somx2*x*l*p[l][l];
  }
}
