#include <math.h>
#include <numerical.h>
#include <vector.h>
#include <genpoly.h>
#include <interp.h>


static double nn, mm;

static void deriv(double xi, double* y, double* f)
{
  if (y[1] < 0.0) {
    f[1] = y[2];
    f[2] = -2.0*y[2]/xi;
  }
  else {
    f[1] = y[2];
    f[2] = -2.0*y[2]/xi - pow(xi,2*mm)*pow(y[1],nn+mm);
  }
}


// Compute generalized polytropes
// normalized so that M=G=1, E_%tot=-1/4

GeneralizedPolytrope::GeneralizedPolytrope(void)
{
  n = -99;
  m = -99;
}


GeneralizedPolytrope::GeneralizedPolytrope(int num, double N, double M, 
					   double eps0, double step)
{

  n = N;
  m = M;
  nn = N;
  mm = M;

  double       by[2], byy[2], bf[2];
  double       *y, *yy, *f, *yl;
  double       eps, a1, a2, a3, xi, fbeg, thetalast, xilast, xistep;
  double       fac, a, b, c, d;
  double       xmassfac, potfac, densfac, rfac;
  int          idim,j,icnt;
      
  y = by-1;
  yy = byy-1;
  f = bf-1;

				/* Compute initial value */
 
  a1=-1.0/(2.0*m + 2.0)/(2.0*m + 3);
  a2=-(n + m)*a1/(4.0*m + 4)/(4.0*m + 5);
  a3=-(n + m)*a2/(6.0*m + 6)/(6.0*m + 7);
 
  xi=step;

  y[1]=1.0 + a1*pow(xi,2*m + 2) + a2*pow(xi,4*m + 4) +
    a3*pow(xi,6*m + 6);

  y[2]=(2.0*m + 2.0)*a1*pow(xi,2*m + 1) +
    (4.0*m + 4.0)*a2*pow(xi,4*m + 3) +
      (6.0*m + 6.0)*a3*pow(xi,6*m + 5);

  fbeg=(2.0*m + 2.0)*(2.0*m + 1.0)*a1*pow(xi,2*m) +
    (4.0*m + 4.0)*(4.0*m + 3.0)*a2*pow(xi,4*m + 2) +
      (6.0*m + 6.0)*(6.0*m + 5.0)*a3*pow(xi,6*m + 4);
 
				/* Begin */
  idim = 2;
  eps = eps0;
  while (y[1] > 0.0) {
    thetalast=y[1];
    deriv(xi, y, yy);
    rk4(y, yy, idim, xi, eps, y, deriv);
    xi += eps;
  }
				/* Iteration */
  xilast=xi-eps;
  icnt=0;
  while( (fabs(xi-xilast) > 1.0e-8 ||
	  fabs(y[1]) > 1.0e-8 ) && icnt < 40 ) {
    xistep=-y[1]*(xilast-xi)/(thetalast-y[1]);
    xilast=xi;
    thetalast=y[1];
    deriv(xi, y, yy);
    rk4(y, yy, idim, xi, xistep, y, deriv);
    xi += xistep;
    icnt++;
  }

  fac=(3.0*m - n + 5.0)/(4.0*m + 6.0);
  densfac=pow(fac,3.0)*7.95774715459476678844e-2* xi/(-y[2]);
  rfac=1.0/(fac*xi);
  xmassfac=-1.0/(xi*xi*y[2]);
  potfac = 1.0/(xi*y[2]);

  eps = (xi-step)/(num-1);
  int nsub = (int)(eps/eps0);
  nsub = nsub>0 ? nsub : 1;
  eps /= nsub;

  dens.x.setsize(1, num);
  dens.y.setsize(1, num);
  dens.y2.setsize(1, num);

  mass.x.setsize(1, num);
  mass.y.setsize(1, num);
  mass.y2.setsize(1, num);

  pot.x.setsize(1, num);
  pot.y.setsize(1, num);
  pot.y2.setsize(1, num);

  dens.x[num] = mass.x[num] = pot.x[num] = rfac*xi;
  dens.y[num] = 0.0;
  mass.y[num] = -xmassfac * xi*xi*y[2];
  pot.y[num] = -fac;

  dens.x[1] = mass.x[1] = pot.x[1] = 0.0;
  if (m==0) 
    dens.y[1] = densfac;
  else
    dens.y[1] = 0.0;
  mass.y[1] = 0.0;
  pot.y[1] = -fac*(1.0-potfac);

  xi = step;

  y[1]=1.0 + a1*pow(xi,2*m + 2) + a2*pow(xi,4*m + 4) +
    a3*pow(xi,6*m + 6);

  y[2]=(2.0*m + 2.0)*a1*pow(xi,2*m + 1) +
    (4.0*m + 4.0)*a2*pow(xi,4*m + 3) +
      (6.0*m + 6.0)*a3*pow(xi,6*m + 5);

  fbeg=(2.0*m + 2.0)*(2.0*m + 1.0)*a1*pow(xi,2*m) +
    (4.0*m + 4.0)*(4.0*m + 3.0)*a2*pow(xi,4*m + 2) +
      (6.0*m + 6.0)*(6.0*m + 5.0)*a3*pow(xi,6*m + 4);
 
  for (int i=2; i<num; i++) {

    for (j=1; j<=nsub; j++) {
      deriv(xi, y, yy);
      rk4(y, yy, idim, xi, eps, y, deriv);
      xi += eps;
    }
    
    dens.x[i] = mass.x[i] = pot.x[i] = rfac*xi;

    dens.y[i] = densfac * pow(xi, 2*m) * pow(y[1],n+m); 

    mass.y[i] = -xmassfac * xi*xi*y[2];

    pot.y[i] = -fac*(1.0 - y[1]*potfac);

  }

  Spline(dens.x, dens.y, 0.0, -1.0e30, dens.y2);
  Spline(mass.x, mass.y, 0.0, -1.0e30, mass.y2);
  Spline(pot.x, pot.y, 0.0, -1.0e30, pot.y2);
  dens.num = mass.num = pot.num = num;

  KF = densfac/pow(-fac*potfac, n+m)/pow(rfac, 2.0*m)
    * 0.5/M_PI / pow(2.0,m-0.5) *
    exp(lgamma(1.0+n+m) + lgamma(0.5+m+n) - lgamma(0.5) - lgamma(0.5+n+m) -
	lgamma(m+1.0) - lgamma(n-0.5));

  dist_defined = TRUE;
}

double GeneralizedPolytrope::get_mass(double r)
{
  double ret;
  if (r>mass.x[mass.num]) return mass.y[mass.num];
  Splint1(mass.x, mass.y, mass.y2, r, ret, 1);
  return ret;
}

double GeneralizedPolytrope::get_density(double r)
{
  double ret;
  if (r>dens.x[dens.num]) return 0.0;
  Splint1(dens.x, dens.y, dens.y2, r, ret, 1);
  return ret;
}

double GeneralizedPolytrope::get_pot(double r)
{
  double ret;
  if (r>dens.x[dens.num]) return -mass.y[mass.num]/r;
  Splint1(pot.x, pot.y, pot.y2, r, ret, 1);
  return ret;
}

double GeneralizedPolytrope::get_dpot(double r)
{
  double ret, dum;
  if (r>dens.x[dens.num]) return mass.y[mass.num]/(r*r);
  Splint2(pot.x, pot.y, pot.y2, r, dum, ret, 1);
  return ret;
}

void GeneralizedPolytrope::get_pot_dpot(double r, double& potl, double& dpotl)
{
  if (r>dens.x[dens.num]) {
    potl = -mass.y[mass.num]/r;
    dpotl = mass.y[mass.num]/(r*r);
  }
  else
    Splint2(pot.x, pot.y, pot.y2, r, potl, dpotl, 1);
}

double GeneralizedPolytrope::get_dpot2(double r)
{
  double ret, den, dpot, dum;
  Splint1(dens.x, dens.y, dens.y2, r, den, 1);
  Splint2(pot.x, pot.y, pot.y2, r, dum, dpot, 1);

  if (r<=0.0)
    ret = 4.0*M_PI*den;
  else
    ret = 4.0*M_PI*den - 2.0*dpot/r;

  return ret;
}

double GeneralizedPolytrope::distf(double E, double L)
{
  if (E>pot.y[pot.num]) return 0.0;
  return KF * pow(pot.y[pot.num] - E, n-1.5) * pow(L, 2.0*m);
}

double GeneralizedPolytrope::dfde(double E, double L)
{
  if (E>pot.y[pot.num]) return 0.0;
  return (1.5 - n) * KF * pow(pot.y[pot.num] - E, n-2.5) * pow(L, 2.0*m);
}

double GeneralizedPolytrope::dfdl(double E, double L)
{
  if (E>pot.y[pot.num]) return 0.0;
  return 2.0*m * KF * pow(pot.y[pot.num] - E, n-2.5) * pow(L, 2.0*m-1.0);
}

double GeneralizedPolytrope::d2fde2(double E, double L)
{
  if (E>pot.y[pot.num]) return 0.0;
  return (1.5 - n) * (2.5 - n) * KF * pow(pot.y[pot.num] - E, n-3.5) * 
    pow(L, 2.0*m);
}

