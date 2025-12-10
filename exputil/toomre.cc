//
// Toomre's Model 1
// with Kalnajs' m distribution function derived from inverting the
// defining integral equation with a separable form: e^m g(x)
//

#include <cmath>
#include <cstdlib>

#include "massmodel.H"
#include "toomre.H"

const double TOL=1.0e-10;
const double TOL2=1.0e-14;
const int ITMAX=20000;

void ToomreDisk::pdist(double E, double L)
{
  static double ee = 1.0;
  static double ll = -1.0;

  if (fabs(e-ee) < TOL && fabs(L-ll) < TOL) return;

  ee = E;
  ll = L;

  e = - E;
  x = sqrt(2*e)*L;

  double cur=1.0, fac;
  double logx = log(x);

  p0 = p1 = p2 = 0.0;

  int j;

  for (j=0; j<ITMAX; j++) {

    if (j>0 && fabs(cur/p0) < TOL2) break;

    fac = lgamma(0.5*(1 + m) + j) - lgamma(0.5*(1 + m)) +
          lgamma(0.5*m + 1.0 + j) - lgamma(0.5*m + 1.0) + 
	  lgamma(0.5*m - 1.5 + j) - lgamma(0.5*m - 1.5) -
	  lgamma(0.5 + j)         + lgamma(0.5)         -
	  lgamma(m + j)           + lgamma(m)           -
	  lgamma(1.0 + j)         + lgamma(1.0)         ;

    cur = exp(fac + logx*2*j);

    p0 += cur;
    if (j==0) continue;
    p1 += exp(fac + logx*(2*j-1)) * 2*j;
    p2 += exp(fac * logx*(2*j-2)) * 2*j*(2*j-1);

  }

#ifdef DEBUG
  if (j>=ITMAX) {
    cerr << "No convergence for [" << E << "," << L << "]\n";
  }
#endif

}

double ToomreDisk::distf(double E, double L)
{

  pdist(E, L);

  return p0 * pow(e, m-1.0) * m / (4.0*M_PI*M_PI);
  
}

double ToomreDisk::dfde(double E, double L)
{

  pdist(E, L);

  return -(
	  p0 * pow(e, m-2.0) * m*(m-1) +
	  p1 * pow(e, m-1.0) * m * L/sqrt(2.0*e)
	  )/ (4.0*M_PI*M_PI);
  
}


double ToomreDisk::dfdl(double E, double L)
{
  pdist(E, L);

  return p1 * pow(e, m-1.0) * m * sqrt(2.0*e) / (4.0*M_PI*M_PI);
}


double ToomreDisk::d2fde2(double E, double L)
{

  pdist(E, L);

  return (
	  p0 * pow(e, m-3.0) * m*(m-1.0)*(m-2.0) -
	  p1 * pow(e, m-2.0) * m*(m-1.0)*2.0/sqrt(2.0*e) -
	  p2 * pow(e, m-1.0) * m/pow(2.0*e, 1.5) +
	  p1 * pow(e, m-1.0) * m/(2.0*e)
	  )/ (4.0*M_PI*M_PI);
  
}



