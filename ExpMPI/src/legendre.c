#define _BSD_SOURCE
#include <stdio.h>
#include <math.h>
  
void legendre_R(int lmax, double x, double **p)
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
	fprintf(stderr, "legendre_R: p[%d][%d]: pll=%e\n", m, m,
		pll);
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p[m][m];
    p[m+1][m] = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p[l][m] = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (isnan(p[l][m]))
	fprintf(stderr, "legendre_R: p[%d][%d]: pl1=%e pl2=%e\n", l, m,
		pl1, pl2);

      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (isnan(x))
    fprintf(stderr, "legendre_R: x\n");
  for(l=0; l<=lmax; l++)
    for (m=0; m<=lmax; m++)
      if (isnan(p[l][m]))
	fprintf(stderr, "legendre_R: p[%d][%d]\n", l, m);


}


void legendre(int lmax, double x, double ***pp)
{
  double fact, somx2, pll, pl1, pl2, **p;
  int m, l;

  p = *pp;
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

  if (isnan(x))
    fprintf(stderr, "legendre: x\n");
  for(l=0; l<=lmax; l++)
    for (m=0; m<=lmax; m++)
      if (isnan(p[l][m]))
	fprintf(stderr, "legendre: p[%d][%d]\n", l, m);


}


void dlegendre_R(int lmax, double x, double **p, double **dp)
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


void dlegendre(int lmax, double x, double ***pp, double ***dpp)
{
  double fact, somx2, pll, pl1, pl2, **p, **dp;
  int m, l;

  p = *pp;
  dp = *dpp;

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
