// This may look like C code, but it is really -*- C++ -*-

//  Routines for computing biorthonormal pairs based on 
//  Clutton-Brock's ultraspherical series
//
//  MDW 11/13/91 [based on "findwake" by MDW: 12/26/1987]

#include <string>
#include <iostream>
#include <cmath>

#include <biorth.H>

double ultra(int n, double l, double x);

CBSphere::CBSphere(void) : AxiSymBiorth(3) {
  
  BiorthID = "CBSphere";
  rbmin = -1.0;
  rbmax = 1.0;
}


double CBSphere::potl(int nn, int l, double x)
{
  int n;

  if (fabs(x) >= 1.0) return 0.0;

  n = nn-1;
  return pow(1.0 - x*x,0.5*(double)l) * sqrt(1.0 - x) * ultra(n,l,x)/
    pow(2.0,0.5+(double)l);
}


void CBSphere::potl(int nn, int l, double x, Eigen::VectorXd& t)
{
  if (fabs(x) >= 1.0) {
    t *= 0.0;
    return;
  }

  int n = nn-1;
  double pfac = pow(1.0 - x*x,0.5*(double)l) * sqrt(1.0 - x) /
    pow(2.0,0.5+(double)l);


  double a,b,u1,u2,u;
  int j;

  u2 = 1.0;
  u  = 2.0*x*(l+1.0);

  t[1] = pfac * u2;
  if (n==0) return;

  t[2] = pfac * u;
  if (n==1) return;

  for (j=2; j<=n; j++) {
    u1 = u2;
    u2 = u;    
    a = 2.0*x*(l+(double)j)/(double)(j);
    b = -(double)(2.0*l + (double)j)/(double)(j);

    u = a*u2 + b*u1;

    t[j+1] = pfac * u;
  }

  return;

}

double CBSphere::dens(int nn, int l, double x)
{
  int n;

  if (fabs(x) >= 1.0) return 0.0;

  n = nn-1;
  return krnl(n,l)/pow(2.0,2.5+(double)l) *
    pow(1.0 - x*x,0.5*(double)l) * pow(1.0 - x,2.5) * ultra(n,l,x);
}

void CBSphere::dens(int nn, int l, double x, Eigen::VectorXd& t)
{
  if (fabs(x) >= 1.0) {
    t *= 0.0;
    return;
  }

  int n = nn-1;
  double pfac = pow(1.0 - x*x,0.5*(double)l) * pow(1.0 - x, 2.5) 
    / pow(2.0,2.5+(double)l);

  double a,b,u1,u2,u;
  int j;

  u2 = 1.0;
  u  = 2.0*x*(l+1.0);

  t[1] = krnl(0, l) * pfac * u2;
  if (n==0) return;

  t[2] = krnl(1, l) * pfac * u;
  if (n==1) return;

  for (j=2; j<=n; j++) {
    u1 = u2;
    u2 = u;    
    a = 2.0*x*(l+(double)j)/(double)(j);
    b = -(double)(2.0*l + (double)j)/(double)(j);

    u = a*u2 + b*u1;

    t[j+1] = krnl(j, l) * pfac * u;
  }

  return;

}


double CBSphere::krnl(int n, int l)
{
  return 4.0*n*(n+2*l+2) + (2*l+1)*(2*l+3);
}

double CBSphere::norm(int n, int l)
{
  return M_PI * krnl(n,l) * exp( 
		    -log(2.0)*((double)(4*l+4))
		    - lgamma((double)(1+n)) - 2.0*lgamma((double)(1+l) )
		    + lgamma((double)(2*l+n+2))
		    )/(double)(l+n+1);
}


/*------------------------------------------------------------------------
 *                                                                       *
 *      Convert between reduced coordinate                               *
 *                                                                       *
 *                 r^2-1                                                 *
 *          rb =  -------                                                *
 *                 r^2+1                                                 *
 *                                                                       *
 *      and its inverse:                                                 *
 *                                                                       *
 *              (1+rb)^(1/2)                                             *
 *          r = ------------                                             *
 *              (1-rb)^(1/2)                                             *
 *                                                                       *
 *-----------------------------------------------------------------------*/


#define BIG 1.0e30
double CBSphere::rb_to_r(double rb)
{
  if (rb>=1.0) 
    return BIG;
  else
    return sqrt( (1.0+rb)/(1.0-rb) );
}

double CBSphere::d_r_to_rb(double r)
{
  double fac;

  fac = r*r + 1.0;;
  return 4.0*r/(fac*fac);
}

double CBSphere::r_to_rb(double r)
{
  return (r*r-1.0)/(r*r+1.0);
}
#undef BIG

  

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


HQSphere::HQSphere(void) : AxiSymBiorth(3) {
  
  BiorthID = "HQSphere";
  rbmin = -1.0;
  rbmax = 1.0;
}


double HQSphere::potl(int nn, int l, double x)
{
  int n;

  if (fabs(x) >= 1.0) return 0.0;

  n = nn-1;
  return pow(1.0 - x*x,(double)l) * (1.0 - x)* ultra(n,2.0*l+0.5,x)/
    pow(2.0,2.0*l + 1.0);
}


void HQSphere::potl(int nn, int ll, double x, Eigen::VectorXd& t)
{
  if (fabs(x) >= 1.0) {
    t *= 0.0;
    return;
  }

  int n = nn-1;
  double pfac = pow(1.0 - x*x,(double)ll) * (1.0 - x) /
    pow(2.0,2.0*ll + 1.0);

  double l = 2.0*ll + 0.5;

  double a,b,u1,u2,u;
  int j;

  u2 = 1.0;
  u  = 2.0*x*(l+1.0);

  t[1] = pfac * u2;
  if (n==0) return;

  t[2] = pfac * u;
  if (n==1) return;

  for (j=2; j<=n; j++) {
    u1 = u2;
    u2 = u;    
    a = 2.0*x*(l+(double)j)/(double)(j);
    b = -(double)(2.0*l + (double)j)/(double)(j);

    u = a*u2 + b*u1;

    t[j+1] = pfac * u;
  }

  return;
}



double HQSphere::dens(int n, int l, double x)
{
  if (fabs(x) >= 1.0) return 0.0;

  return krnl(n,l)/pow(2.0,2.0*l + 2.0) *
    pow(1.0 - x*x,(double)l-1.0) * pow(1.0 - x,5.0) * ultra(n,2.0*l+0.5,x);
}

void HQSphere::dens(int nn, int ll, double x, Eigen::VectorXd& t)
{
  if (fabs(x) >= 1.0) {
    t *= 0.0;
    return;
  }

  int n = nn-1;
  double pfac = pow(1.0 - x*x,(double)ll-1.0) * pow(1.0 - x,5.0) /
    pow(2.0,2.0*ll + 2.0);

  double l = 2.0*ll + 0.5;

  double a,b,u1,u2,u;
  int j;

  u2 = 1.0;
  u  = 2.0*x*(l+1.0);

  t[1] = krnl(0, ll) * pfac * u2;
  if (n==0) return;

  t[2] = krnl(1, ll) * pfac * u;
  if (n==1) return;

  for (j=2; j<=n; j++) {
    u1 = u2;
    u2 = u;    
    a = 2.0*x*(l+(double)j)/(double)(j);
    b = -(double)(2.0*l + (double)j)/(double)(j);

    u = a*u2 + b*u1;

    t[j+1] = krnl(j, ll) * pfac * u;
  }

  return;
}


double HQSphere::krnl(int n, int l)
{
  return 0.5*n*(n+4*l+3) + (l+1)*(2*l+1);
}

double HQSphere::norm(int n, int l)
{
  return M_PI * krnl(n,l) * exp( 
	    -log(2.0)*((double)(8*l+4))
	    - lgamma((double)(1+n)) - 2.0*lgamma((double)(1.5+2.0*l))
		    + lgamma((double)(4*l+n+3))
		    )/(double)(2*l+n+1.5);

}

/*------------------------------------------------------------------------
 *                                                                       *
 *      Convert between reduced coordinate                               *
 *                                                                       *
 *                 r - 1                                                 *
 *          rh =  -------                                                *
 *                 r + 1                                                 *
 *                                                                       *
 *      and its inverse:                                                 *
 *                                                                       *
 *              (1+rh)
 *          r = ------
 *              (1-rh)
 *                                                                       *
 *-----------------------------------------------------------------------*/


#define BIG 1.0e30
double HQSphere::rb_to_r(double rb)
{
  if (rb>=1.0) 
    return BIG;
  else
    return (1.0+rb)/(1.0-rb);
}

double HQSphere::d_r_to_rb(double r)
{
  double fac;

  fac = r + 1.0;
  return 2.0/(fac*fac);
}

double HQSphere::r_to_rb(double r)
{
  return (r-1.0)/(r+1.0);
}
#undef BIG



//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


Eigen::VectorXd sbessjz(int, int);
double jn_sph(int, double);

BSSphere::BSSphere(double RMAX, int NMAX, int LMAX) : AxiSymBiorth(3) {
  
  BiorthID = "BSSphere";

  
  rmax = RMAX;
  nmax = NMAX;
  lmax = LMAX;

  a.resize(lmax+1);

  for (int l=0; l<=lmax; l++) {
    a[l] = sbessjz(l-1,nmax);
  }

  t_n = 0;

}

BSSphere::~BSSphere(void)
{
  // Nothing
}


double BSSphere::potl(int n, int l, double r)
{
  if (l > lmax) bomb("potl: l too large");
  if (n > nmax) bomb("potl: n too large");

  return M_SQRT2/fabs(a[l][n]*jn_sph(l,a[l][n])) * pow(rmax,-0.5) *
    jn_sph(l, a[l][n]*r/rmax);
}


void BSSphere::setup_potl_table(void)
{
  t_n = 40000;
  t_dr.resize(lmax+1);
  t_f.resize(lmax+1, nmax);
  t_g.resize(lmax+1, nmax);
  t_y.resize(lmax+1, t_n+1);

  for (int l=0; l<=lmax; l++) {

    t_dr[l] = a[l][nmax]/t_n;

    int i;
    for (i=0; i<=t_n; i++) 
      t_y(l, i) = jn_sph(l, t_dr[l]*i);

    for (i=1; i<=nmax; i++) {
      t_f(l, i) = M_SQRT2/fabs(a[l][i]*jn_sph(l,a[l][i])) * pow(rmax,-0.5);
      t_g(l, i) = M_SQRT2*fabs(a[l][i]/jn_sph(l,a[l][i])) * pow(rmax,-2.5);
    }
  }

}

void BSSphere::potl(int n, int l, double r, Eigen::VectorXd& t)
{
  if (l > lmax) bomb("potl: l too large");
  if (n > nmax) bomb("potl: n too large");

  if (!t_n) setup_potl_table();

  double r0, rs;
  int indx;
 
  for (int i=1; i<=nmax; i++) {
    rs = a[l][i]*r/rmax;
    indx = (int)(rs/t_dr[l]);
    indx = indx>=t_n ? t_n-1 : indx;

    r0 = t_dr[l]*indx;
   
    t[i] = t_f(l, i)*(t_y(l, indx)*(r0 + t_dr[l] - rs) + 
		      t_y(l, indx+1)*(rs - r0))/t_dr[l];
  }

}


double BSSphere::dens(int n, int l, double r)
{
  if (l > lmax) bomb("dens: l too large");
  if (n > nmax) bomb("dens: n too large");

  return a[l][n]*M_SQRT2/fabs(jn_sph(l,a[l][n])) * pow(rmax,-2.5) * 
    jn_sph(l, a[l][n]*r/rmax);

}

void BSSphere::dens(int n, int l, double r, Eigen::VectorXd& t)
{
  if (l > lmax) bomb("dens: l too large");
  if (n > nmax) bomb("dens: n too large");

  if (!t_n) setup_potl_table();

  double r0, rs;
  int indx;
 
  for (int i=1; i<=nmax; i++) {
    rs = a[l][i]*r/rmax;
    indx = (int)(rs/t_dr[l]);
    indx = indx>=t_n ? t_n-1 : indx;

    r0 = t_dr[l]*indx;
   
    t[i] = t_g(l, i)*(t_y(l, indx)*(r0 + t_dr[l] - rs) + 
		      t_y(l, indx+1)*(rs - r0))/t_dr[l];
  }

}


