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

#include <functional>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "OrthoPoly.H"
#include "CauchyPV.H"

double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void tqli(Eigen::VectorXd& d,
	  Eigen::VectorXd& e,
	  std::vector<Eigen::VectorXd>& z, int n)
{
  double pythag(double a, double b);

  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  for (int i=1; i<=n; i++) e[i-1] = e[i];
  e[n] = 0.0;

  for (int l=1; l<=n; l++) {
    iter=0;
    do {
      for (int m=l; m<=n-1; m++) {
	dd = fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  std::cerr << "Too many iterations in tqli" << std::endl;
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

				// Constructor
PVQuad :: PVQuad(int N)
{
  n = N;

  chrs.resize(n+1);
  
  get_eigenvalues();
  get_orthopolies();
}

void PVQuad::get_eigenvalues(void)
{
  std::vector<Eigen::VectorXd> z(n+1);

  for (int i=0; i<=n; i++) {
    z[i].resize(n+1);
    z[i].setZero();
    z[i][i] = 1.0;
  }

  Eigen::VectorXd work1(n+1);
  Eigen::VectorXd work2(n+1);
  
  work1[0] = -coef2(0)/coef3(0);

  for (int i=0; i<n; i++) {
    work1[i+1] = -coef2(i)/coef3(i);
    work2[i+1] = sqrt(coef4(i)*coef1(i-1)/(coef3(i-1)*coef3(i)));
  }

  tqli(work1, work2, z, n+1);

  eign.resize(n+1);
  chrs.resize(n+1);
 
  for (int i=0; i<=n; i++) eign[i] = work1[i+1];
  for (int i=0; i<=n; i++)
    chrs[i] = z[i+1][1]*z[i+1][1]/(0.5*z[i+1].adjoint()*z[i+1]);

}

void PVQuad::get_orthopolies(void)
{
  Eigen::VectorXd temp(n+1);
  poly.resize(n+1);

  for (auto & v : poly) v.resize(n+1);

  for (int i=0; i<=n; i++) {
    temp = fv(eign[i], n);
    for (int j=0; j<=n; j++) poly[j][i] = temp[j];
  }

				// Normalize
  nrml.resize(n+1);
  Eigen::VectorXd comb(n+1);
  for (int j=0; j<=n; j++) {
    for (int k=0; k<=n; k++) comb[k] = chrs[k] * poly[j][k];
    nrml[j] = comb.adjoint()*poly[j];
  }
}

Eigen::VectorXd& PVQuad::return_coefs(std::function<double(double)> func)
{
  coefs.resize(n+1);

  Eigen::VectorXd input(n+1);
  for (int j=0; j<=n; j++) input[j] = func(eign[j]);

  return get_coefs(input);
}


Eigen::VectorXd& PVQuad::get_coefs(Eigen::VectorXd& input)
{
  coefs.resize(n+1);

  Eigen::VectorXd comb(n+1);
  for (int i=0; i<=n; i++) {
    for (int j=0; j<=n; j++) comb[j] = chrs[j] * poly[i][j];
    coefs[i] = comb.dot(input) / nrml[i];
  }

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





