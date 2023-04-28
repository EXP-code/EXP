#include <functional>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <cmath>

#include <Bess.H>
#include <EXPmath.H>

/**
   Brent's algorithm root finder
*/  
double ZBrent(std::function<double(double)> f, double x1, double x2, double tol)
{
  const int    ITMAX = 500;
  const double EPS   = 1.0E-15;

  double a=x1, b=x2, c, d, e;

  double fa=f(a), fb=f(b), fc;

  if (fb*fa > 0.0) {
    throw std::runtime_error("Root must be bracketed in ZBRENT");
  }
  fc = fb;

  for (int iter=0; iter<ITMAX; iter++) {

    if (fb*fc > 0.0) {
      c  = a;
      fc = fa;
      e  = d = b-a;
    }

    if (fabs(fc) < fabs(fb)) {
      a  = b;
      b  = c;
      c  = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    double tol1 = 2.0*EPS*fabs(b)+0.5*tol;
    double xm = 0.5*(c-b);

    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e)  >= tol1 && fabs(fa) > fabs(fb)) {
      double s = fb/fa, p, q;
      if (a == c) {
	p = 2.0*xm*s;
	q = 1.0-s;
      } else {
	q = fa/fc;
	double r = fb/fc;
	p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p = fabs(p);
      double min1 = 3.0*xm*q-fabs(tol1*q);
      double min2 = fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e = d;
	d = p/q;
      } else {
	d = xm;
	e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a  = b;
    fa = fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    fb = f(b);
  }
  std::cerr << "Maximum number of iterations exceeded in ZBRENT";

  return b;
}

/**
   Find bracket for prior to refinement
*/  
std::tuple<double, double>
Bess::bracket(std::function<double(double)> B, double zmin, double z0)
{
  double del = 0.1;
  double cur = B(z0);
  double  zl = z0;
  double   z = zl - del;
  
  bool found = false;
  while (z > zmin) {
    double tst = B(z);
    if (tst*cur <= 0.0) {
      found = true;
      break;
    }
    zl = z;
    z -= del;
    cur = tst;
  }
    
  if (not found) {

    cur = B(z0);
    zl  = z0;
    z   = zl + del;
    
    while (z < z0+0.5*M_PI) {
      double tst = B(z);
      if (tst*cur <= 0.0) {
	found = true;
	break;
      }
      zl = z;
      z += del;
      cur = tst;
    }
  }

  if (not found) throw std::runtime_error("Could not find valid bracket");

  if (zl<=z) return {zl, z};
  return {z, zl};
}

std::vector<double> Bess::zeros(double Nu, int Nz, double tol)
{
  std::function<double(double)> S = [](double x){return x*x;  };
  std::function<double(double)> C = [](double x){return x*x*x;};

  std::vector<double> approx(Nz);

  // Compute asymptotic root for starting point
  //
  for(int k=0; k<Nz; k++) {
    double mu = 4.0*S(Nu);
    double beta = (static_cast<double>(k+1)+0.5*Nu-0.25)*M_PI;
    approx[k]  = beta - (mu-1.0)/8./beta - 4.0*(mu-1.0)*(7.0*mu-31.0)/3.0/C(8.0*beta);
    approx[k] -=  32.*(mu-1.0)*(83.0*S(mu)-982*mu+3779.0)/15./C(8.0*beta)/S(8.0*beta);
    approx[k] -= 64.*(mu-1.0)*(6949.*C(mu)-153855.0*S(mu)+1585743*mu-6277237)/105./S(C(8.0*beta))/(8.0*beta); 

  }

  if (log10(1.0/tol/approx.back())<2.0) {
    std::cerr << "Warning: you want more roots than accuracy can guarantee" << std::endl;
  }
    

  std::vector<double> refine(Nz);

  std::function<double(double)> B = [Nu](double x){return EXPmath::cyl_bessel_j(Nu, x);};

  for (int k=0; k<Nz; k++) {
    try {
      std::tuple<double, double> R;
      if (k==0)
	R = bracket(B, 0.0, approx[k]);
      else
	R = bracket(B, refine[k-1], approx[k]);

      refine[k] = ZBrent(B, std::get<0>(R), std::get<1>(R), tol);
    }
    catch (std::runtime_error& e) {
      refine[k] = approx[k];
    }
  }

  return refine;
}
