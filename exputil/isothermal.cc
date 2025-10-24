// This may look like C code, but it is really -*- C++ -*-

#include <cstdlib>
#include <limits>
#include <cmath>

#include <iostream>
#include <fstream>

#include <string>
#include <vector>

#include "numerical.H"
#include "interp.H"
#include "isothermal.H"

static double sig;
static void iso_derivs(double x, Eigen::VectorXd& y, Eigen::VectorXd& dy)
{
  dy[1] = y[2];
  double den = 9.0*sig*sig*exp(-y[1]/(sig*sig));
  if (x<=0.0)
    dy[2] = 3.0*sig*sig;
  else
    dy[2] = den - 2.0*y[2]/x;
  dy[3] = den*x*x;
}


IsothermalSphere::IsothermalSphere(double RCORE, double RMAX, double VROT,
				   int NUM, int dN)
{
  vrot = VROT;
  sigma = 1.0;
  rscale = RCORE; 
  rmax = RMAX;
  dist_defined = true;
  ModelID = "IsothermalSphere"; 
  dim = 3;

  r.resize(NUM);
  d.resize(NUM);
  m.resize(NUM);
  p.resize(NUM);
  d2.resize(NUM);
  m2.resize(NUM);
  p2.resize(NUM);

				// Begin integration

  double tol = 1.0e-10;
  double dr = RMAX/RCORE/NUM;
  double h = dr/dN;

  Eigen::VectorXd y(3);

				// Initial conditions
  r[0] = 0.0;
  d[0] = 9.0*sigma*sigma/(4.0*M_PI);
  m[0] = 0.0;
  p[0] = 0.0;

  F = d[1]/(4.0*M_PI*sqrt(M_PI/2.0)/sigma*sigma*sigma);

  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;

  sig = sigma;

  for (int i=1; i<NUM; i++) {
    
    r[i] = r[i-1] + dr;

    integrate_ode(y, r[i-1], r[i], h, tol, 3, iso_derivs, rkqc);

    d[i] = d[0] * exp(-y[0]/(sigma*sigma));
    p[i] = y[0];
    m[i] = y[2];
  }

				// Scaling: vrot at edge

  double beta = r[NUM]/m[NUM] * vrot*vrot;
  double gamma = 1.0/(rscale*rscale*sqrt(beta));
  double dfac = gamma * pow(beta, 1.5);

  r *= rscale;
  d *= dfac;
  m *= rscale*rscale*rscale*dfac;
  p *= beta;

				// Set up splines

  Spline(r, d, 0.0, -1.0e30, d2);
  Spline(r, m, 0.0, -1.0e30, m2);
  Spline(r, p, 0.0, -1.0e30, p2);
}


//======================================================================


static double fconst;
static double pmax;

static void sl_derivs(double x, Eigen::VectorXd& y, Eigen::VectorXd& dy)
{
  dy[1] = y[2];

  double vm = sqrt(2.0*fabs(pmax - y[1]));
  double norm1 = 4.0*M_PI*
    ( - vm*exp(-0.5*vm*vm) + sqrt(0.5*M_PI) * erf(vm/sqrt(2.0)) );
  double norm2 = 4.0*M_PI*vm*vm*vm/3.0;

  double dm = (norm1*exp(-y[1]) - norm2*fconst)*exp(2*x);
  dy[2] = dm - y[2];
  dy[3] = dm * exp(x);
  dy[4] = 0.5*dm*y[1] * exp(x);
}


LowSingIsothermalSphere::LowSingIsothermalSphere
(double RMIN, double RMAX, int NUM, int dN)
{
  num = NUM;
  dist_defined = true;
  ModelID = "LowSingIsothermal"; 
  dim = 3;

  u.resize(NUM);
  d.resize(NUM);
  m.resize(NUM);
  p.resize(NUM);
  d2.resize(NUM);
  m2.resize(NUM);
  p2.resize(NUM);

				// Begin integration
  double tol = 1.0e-10;
  double umin = log(RMIN);
  double umax = log(RMAX);
  double du = (umax - umin)/NUM;
  double h = du/dN;

  Eigen::VectorXd y(4);

				// Initial conditions
  y[0] = 2.0*umin;
  y[1] = 2.0;
  y[2] = exp(umin);
  y[3] = (umin - 1.0)*exp(umin);

  pmax = 2.0*umax;
  fconst = exp(-pmax);

  double ucur, ulast;

  ucur = umin;
  while (y[0] < pmax) {
    ulast = ucur;
    ucur += du;
    integrate_ode(y, ulast, ucur, h, tol, 4, sl_derivs, rkqc);
  }

  while (fabs(pmax - y[0]) > 1.0e-8) {
    ulast = ucur;
    du = (pmax - y[0])/y[1];
    ucur = ulast + du;
    h = du/dN;
    integrate_ode(y, ulast, ucur, h, tol, 4, sl_derivs, rkqc);
  }

				// Regrid solution
  du = (ucur - umin)/(NUM-1);
  h = du/dN;
  F = 1.0/(4.0*M_PI);

  y[0] = 2.0*umin;
  y[1] = 2.0;
  y[2] = exp(umin);
  y[3] = (umin - 1.0)*exp(umin);

  double vm = sqrt(2.0*fabs(pmax - y[1]));
  double norm1 = 4.0*M_PI*
    ( - vm*exp(-0.5*vm*vm) + sqrt(0.5*M_PI) * erf(vm/sqrt(2.0)) );
  double norm2 = 4.0*M_PI*vm*vm*vm/3.0;
  
  u[0] = umin;
  d[0] = F * (norm1*exp(-y[1]) - norm2*fconst);
  m[0] = exp(umin);
  p[0] = 2.0*umin;

  for (int i=1; i<NUM; i++) {
    
    u[i] = u[i-1] + du;

    integrate_ode(y, u[i-1], u[i], h, tol, 4, sl_derivs, rkqc);

    vm = sqrt(2.0*fabs(pmax - y[0]));
    norm1 = 4.0*M_PI*
    ( - vm*exp(-0.5*vm*vm) + sqrt(0.5*M_PI) * erf(vm/sqrt(2.0)) );
    norm2 = 4.0*M_PI*vm*vm*vm/3.0;

    d[i] = F * (norm1*exp(-y[0]) - norm2*fconst);
    p[i] = y[0];
    m[i] = y[2];
  }


  double pfac = p[num-2] + m[num-2]*exp(-u[num-2]);
  for (int n=0; n<num; n++) p[n] -= pfac;

  double M0 = y[1];
  double W0 = y[2] - 0.5*y[1]*pfac;

  double beta = -0.5*M0/W0;
  double gamma = pow(beta, 3.5) * pow(-2.0*W0, 2.0);
  double dfac = gamma * pow(beta, 1.5);
  double rfac = pow(beta, -0.25) * pow(gamma, -0.5);

  for (int n=0; n<num; n++) {
    u[n] += log(rfac);
    d[n] *= dfac;
    m[n] /= y[2];
    p[n] *= beta;
  }

  sigma = sqrt(beta);
  rmin = exp(u[1]);
  rmax = exp(u[num]);

  F *= gamma * exp(-pfac);
  Fconst = F * exp(-p[num]/(sigma*sigma));

				// Set up splines

  Spline(u, d, 0.0, -1.0e30, d2);
  Spline(u, m, 0.0, -1.0e30, m2);
  Spline(u, p, 0.0, -1.0e30, p2);

}

