#include <cmath>
#include "expand.H"
#include <Hernquist.H>

Hernquist::Hernquist(Component* c0, const YAML::Node& conf, MixtureBasis* m) :
  SphericalBasis(c0, conf, m)
{
  id = "Hernquist sphere";
  initialize();
  setup();
}

void Hernquist::initialize(void)
{
  // Do nothing . . .
}

double Hernquist::knl(int n, int l)
{
  return 0.5*n*(n+4*l+3) + (l+1)*(2*l+1);
}

double Hernquist::norm(int n, int l)
{
  extern double dgammln(double);
  return M_PI * knl(n, l) * 
    exp( 
	-log(2.0)*((double)(8*l+4))
	- lgamma((double)(1+n)) - 2.0*lgamma((double)(1.5+2.0*l))
	+ lgamma((double)(4*l+n+3))
	)/(double)(2*l+n+1.5);
}

void Hernquist::get_dpotl(int lmax, int nmax, double r,
			  Eigen::MatrixXd& p, Eigen::MatrixXd& dp,
			  int tid)
{
  double x     = r_to_xi(r);
  double dx    = d_r_to_xi(r);
  double fac   = 0.25*(1.0 - x*x);
  double rfac  = 0.5*(1.0 - x);
  double drfac = -1.0/(1.0 - x*x);
  
  for (int l=0; l<=lmax; l++) {
    double dfac1 = 1.0 + x + 2.0*x*l;
    double dfac2 = 4.0*l + 3.0;

    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    get_ultra(nmax-1, 2.0*l+1.5, x, du[tid]);

    for (int n=0; n<nmax; n++) p(l, n) = rfac*u[tid][n];
    dp(l, 0) = dx*drfac*dfac1*rfac*u[tid][0];
    for (int n=1; n<nmax; n++)
      dp(l, n) = dx*rfac*(drfac*dfac1*u[tid][n] + dfac2*du[tid][n-1]);

    rfac *= fac;
  }

}

void Hernquist::get_potl(int lmax, int nmax, double r,
			 Eigen::MatrixXd& p, int tid)
{
  double x    = r_to_xi(r);
  double fac  = 0.25*(1.0 - x*x);
  double rfac = 0.5*(1.0 - x);
  
  for (int l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    for (int n=0; n<nmax; n++) p(l, n) = rfac*u[tid][n];
    rfac *= fac;
  }
}

void Hernquist::get_dens(int lmax, int nmax, double r,
			 Eigen::MatrixXd& p, int tid)
{ 
  double x    = r_to_xi(r);
  double fac  = 0.25*(1.0 - x*x);
  double rfac = 0.25*pow(1.0 - x, 5.0)/(1.0 - x*x);

  for (int l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    for (int n=0; n<nmax; n++) p(l, n) = krnl(l, n)*rfac*u[tid][n];
    rfac *= fac;
  }

}

void Hernquist::get_potl_dens(int lmax, int nmax, double r, 
			      Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid)
{
  double x     = r_to_xi(r);
  double fac   = 0.25*(1.0 - x*x);
  double rfacp = 0.5*(1.0 - x);
  double rfacd = 0.25*pow(1.0 - x, 5.0)/(1.0 - x*x);


  for (int l=0; l<=lmax; l++) {
    get_ultra(nmax-1, 2.0*l+0.5, x, u[tid]);
    for (int n=0; n<nmax; n++) p(l, n) = rfacp*u[tid][n];
    for (int n=0; n<nmax; n++) d(l, n) = krnl(l, n)*rfacd*u[tid][n];
    rfacp *= fac;
    rfacd *= fac;
  }
}

/*------------------------------------------------------------------------
 *                                                                       *
 *      Convert between reduced coordinate                               *
 *                                                                       *
 *               r - 1                                                   *
 *          x = -------                                                  *
 *               r + 1                                                   *
 *                                                                       *
 *      and its inverse:                                                 *
 *                                                                       *
 *              (1+x)                                                    *
 *          r = ------                                                   *
 *              (1-x)                                                    *
 *                                                                       *
 *-----------------------------------------------------------------------*/


#define BIG 1.0e30
double Hernquist::xi_to_r(double x)
{
  if (x>=1.0) 
    return BIG;
  else
    return (1.0+x)/(1.0-x);
}

double Hernquist::d_r_to_xi(double r)
{
  double fac;

  fac = r + 1.0;
  return 2.0/(fac*fac);
}

double Hernquist::r_to_xi(double r)
{
  return (r-1.0)/(r+1.0);
}
#undef BIG
