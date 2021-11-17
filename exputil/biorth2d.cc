//  Routines for computing biorthonormal pairs based on 
//  Clutton-Brock's 2-dimensional series
//

#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>

#include <biorth2d.H>
#include <OrthoPoly.H>
#include <gaussQ.H>


CBDisk::CBDisk(void) : AxiSymBiorth(2) {
  BiorthID = "CBDisk";
  numz = 80;
}


double CBDisk::phif(int np, int m, double r)
{
  int n = np-1;
				// By recurrance relation

  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double cur = sqrt(fac);

  for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);

  if (n==0) return cur;
  
  double curl1 = cur;
  double curl2;
  
  fac *= r2 - 1.0;
  cur *= fac*(2*m+1);

  if (n==1) return cur;

  for (int nn=2; nn<=n; nn++) {
    curl2 = curl1;
    curl1 = cur;
    cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
      (1.0 + (double)(2*m-1)/nn)*curl2;
  }

  return cur;
}

double CBDisk::potl(int np, int m, double r)
{
  return pow(r, (double)m+1.0e-20) * phif(np, m, r);
}


				// By recurrance relation

void CBDisk::potl(int np, int m, double r, Eigen::VectorXd& a)
{
  int n = np-1;
  double pfac = pow(r, (double)m+1.0e-20);
  
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double cur = sqrt(fac);

  for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);

  a[1] = pfac*cur;

  if (n==0) return;
  
  double curl1 = cur;
  double curl2;
  
  fac *= r2 - 1.0;
  cur *= fac*(2*m+1);

  a[2] = pfac*cur;

  if (n==1) return;

  for (int nn=2; nn<=n; nn++) {
    curl2 = curl1;
    curl1 = cur;
    cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
      (1.0 + (double)(2*m-1)/nn)*curl2;
    a[nn+1] = pfac*cur;
  }
  
  return;
}



double CBDisk::dens(int nn, int m, double r)
{

  if (nn>2) 
    return pow(r, (double)m+1.0e-20)*(phif(nn, m+1, r) - phif(nn-2, m+1, r));
  else
    return pow(r, (double)m+1.0e-20)*phif(nn, m+1, r);
}


void CBDisk::dens(int np, int mm, double r, Eigen::VectorXd& a)
{
  int n = np-1;
  double pfac = pow(r, (double)mm+1.0e-20);
  
  int m = mm+1;
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double cur = sqrt(fac);

  for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);

  a[1] = pfac*cur;

  if (n==0) return;
  
  double curl1 = cur;
  double curl2;
  
  fac *= r2 - 1.0;
  cur *= fac*(2*m+1);

  a[2] = pfac*cur;

  if (n==1) return;

  int nn;
  for (nn=2; nn<=n; nn++) {
    curl2 = curl1;
    curl1 = cur;
    cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
      (1.0 + (double)(2*m-1)/nn)*curl2;
    a[nn+1] = pfac*cur;
  }
  
  for (nn=n; nn>2; nn--)
    a[nn] -= a[nn-2];
  
  return;
}


double CBDisk::norm(int n, int m)
{
  double ans = 1.0;
  
  for (int i=n+1; i<=n+2*m; i++)
    ans *= i;

  return pow(0.5, 2*m+1)*ans;
}


double CBDisk::potlRZ(int np, int m, double r, double z)
{
  GenLagu lag(2.0*m);
  static int M=-1;
  static std::shared_ptr<LaguQuad> L;
  if (m!=M) {
    L = std::make_shared<LaguQuad>(numz, 2*m);
    M = m;
  };

  double q;
  double ans=0.0;
  double fac=1.0/(1.0+fabs(z));

  for (int i=0; i<numz; i++) {
    q = L->knot(i);
    ans += L->weight(i) * jn(m, fac*r*q) * lag.f(2.0*fac*q, np-1) * 
      pow(q, -m);
  }

//  return ans * fac*pow(fac/r, (double)m);
  return ans * pow(fac, m+1);
}


//  Routines for computing 2d biorthonormal pairs based on 
//  Bessel functions
//

Eigen::VectorXd bessjz(int, int);

BSDisk::BSDisk(double RMAX, int NMAX, int MMAX) : AxiSymBiorth(2)
{
  BiorthID = "BSDisk";
  numz = 80;
  rmax = RMAX;
  nmax = NMAX;
  mmax = MMAX;

  a.resize(mmax+2);

  for (int m=0; m<=mmax+1; m++) {
    a[m] = bessjz(m, nmax);
  }

  t_n = 0;

}

BSDisk::~BSDisk(void)
{
}


double BSDisk::potl(int n, int m, double r)
{
  if (m > mmax) bomb("potl: m too large");
  if (n > nmax) bomb("potl: n too large");

  double nrm = sqrt(1.0 - (double)(m*m)/(a[m][n]*a[m][n])) * jn(m, a[m][n]);

  return M_SQRT2/fabs(rmax*nrm) * jn(m, a[m][n]*r/rmax);
}


void BSDisk::setup_potl_table(void)
{
  t_n = 40000;
  t_dr.resize(mmax);
  t_f. resize(mmax, nmax);
  t_y. resize(mmax, t_n+1);

  double nrm;

  for (int m=0; m<=mmax; m++) {

    t_dr[m] = a[m][nmax-1]/t_n;

    for (int i=0; i<=t_n; i++) 
      t_y(m, i) = jn(m, t_dr[m]*i);

    for (int i=0; i<nmax; i++) {
      nrm = sqrt(1.0 - (double)(m*m)/(a[m][i]*a[m][i])) * jn(m, a[m][i]);
      t_f(m, i) = M_SQRT2/fabs(rmax*nrm);
    }
  }

}

void BSDisk::potl(int n, int m, double r, Eigen::VectorXd& t)
{
  if (m > mmax) bomb("potl: m too large");
  if (n > nmax) bomb("potl: n too large");

  if (!t_n) setup_potl_table();

  double r0, rs;
  int indx;
 
  for (int i=1; i<=nmax; i++) {
    rs = a[m][i]*r/rmax;
    indx = (int)(rs/t_dr[m]);
    indx = indx>=t_n ? t_n-1 : indx;

    r0 = t_dr[m]*indx;
   
    t[i] = t_f(m, i)*(t_y(m, indx)*(r0 + t_dr[m] - rs) + 
		      t_y(m, indx+1)*(rs - r0))/t_dr[m];
  }

}


double BSDisk::dens(int n, int m, double r)
{
  if (m > mmax) bomb("dens: m too large");
  if (n > nmax) bomb("dens: n too large");

  double nrm = sqrt(1.0 - (double)(m*m)/(a[m][n]*a[m][n])) * jn(m, a[m][n]);

  return M_SQRT2/fabs(rmax*nrm) * jn(m, a[m][n]*r/rmax);
}

void BSDisk::dens(int n, int m, double r, Eigen::VectorXd& t)
{
  potl(n, m, r, t);
}

double BSDisk::norm(int n, int m)
{
  return 1.0;
}


double BSDisk::potlRZ(int n, int m, double r, double z)
{
  return potl(n, m, r) * exp(-fabs(z));
}


