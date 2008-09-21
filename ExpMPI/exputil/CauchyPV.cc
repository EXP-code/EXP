/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Compute Cauchy principal value integrals using quadrature rules following
 *  Paget and Elliott (Numer. Math. 19, 373).
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 07/14/94
 *
 ***************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Vector.h>
#include <OrthoPoly.h>
#include <CauchyPV.h>

using namespace std;

				// Constructor
PVQuad :: PVQuad(int N)
{
  n = N;

  chrs.setsize(0, n);
  
  get_eigenvalues();
  get_orthopolies();
}

void PVQuad::get_eigenvalues(void)
{
  void tqli(Vector& a, Vector& b, Vector *z, int n);
  
  Vector* z = new Vector [n+1] - 1;
  int i;

  for (i=1; i<=n+1; i++) {
    z[i].setsize(1, n+1);
    z[i].zero();
    z[i][i] = 1.0;
  }

  Vector work1(1, n+1);
  Vector work2(1, n+1);
  
  work1[1] = -coef2(0)/coef3(0);

  for (i=1; i<=n; i++) {
    work1[i+1] = -coef2(i)/coef3(i);
    work2[i+1] = sqrt(coef4(i)*coef1(i-1)/(coef3(i-1)*coef3(i)));
  }

  tqli(work1, work2, z, n+1);

  eign.setsize(0, n);
  chrs.setsize(0, n);
 
  for (i=0; i<=n; i++) eign[i] = work1[i+1];
  for (i=0; i<=n; i++) chrs[i] = z[i+1][1]*z[i+1][1]/(0.5*z[i+1]*z[i+1]);

}

void PVQuad::get_orthopolies(void)
{
  Vector temp(0, n);
  poly = new Vector [n+1];

  int i;
  for (i=0; i<=n; i++) poly[i].setsize(0, n);

  for (i=0; i<=n; i++) {
    temp = fv(eign[i], n);
    for (int j=0; j<=n; j++) poly[j][i] = temp[j];
  }

				// Normalize
  nrml.setsize(0, n);
  for (int j=0; j<=n; j++) nrml[j] = (chrs & poly[j])*poly[j];
}


Vector& PVQuad::return_coefs(double (*func)(double))
{
  coefs.setsize(0, n);

  Vector input(0, n);
  for (int j=0; j<=n; j++) input[j] = func(eign[j]);

  return get_coefs(input);
}


Vector& PVQuad::get_coefs(Vector& input)
{
  coefs.setsize(0, n);

  for (int i=0; i<=n; i++)
    coefs[i] = ((chrs & poly[i]) * input) / nrml[i];
  
  return coefs;
}

double PVQuad::get_integral(double lambda)
{
  double b0=0.0, b1=0.0, b2;

  for (int k=n; k>=0; k--) {
    b2 = b1;
    b1 = b0;

    b0 = coefs[k] + (coef3(k)/coef1(k)*lambda + coef2(k)/coef1(k)) * b1
      - coef4(k+1)/coef1(k+1) * b2;
  }

  return b0*q0(lambda) - b1*nrml[0]*coef3(0)/coef1(0);
}


//===========================================================================

// From Numerical Recipes


double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


void tqli(Vector& d, Vector& e, Vector *z, int n)
{
  double pythag(double a, double b);
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  for (i=2; i<=n; i++) e[i-1] = e[i];
  e[n] = 0.0;
  for (l = 1; l<=n; l++) {
    iter=0;
    do {
      for (m=l; m<=n-1; m++) {
	dd = fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  cerr << "Too many iterations in tqli\n";
	  exit(-1);
	}
	g = (d[l+1]-d[l])/(2.0*e[l]);
	r = pythag(g,1.0);
	g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s = c = 1.0;
	p = 0.0;
	for (i=m-1;i>=l;i--) {
	  f = s*e[i];
	  b = c*e[i];
	  e[i+1] = (r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d[i+1]-p;
	  r = (d[i]-g)*s+2.0*c*b;
	  d[i+1] = g+(p=s*r);
	  g = c*r-b;
	  for (k=1; k<=n; k++) {
	    f = z[i+1][k];
	    z[i+1][k] = s*z[i][k] + c*f;
	    z[i][k] = c*z[i][k] - s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
}


double pythag(double a, double b)
{
  double absa,absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

