#include "expand.H"
#include <CBrock.H>

CBrock::CBrock(Component* c0, const YAML::Node& conf, MixtureBasis *m) :
  SphericalBasis(c0, conf, m)
{
  id = "Clutton-Brock sphere";
  initialize();
  setup();
}

double CBrock::knl(int n, int l)
{
  return 4.0*n*(n+2*l+2) + (2*l+1)*(2*l+3);
}

double CBrock::norm(int n, int l)
{
  return M_PI * knl(n,l) * 
    exp( 
	-log(2.0)*((double)(4*l+4))
	- lgamma((double)(1+n)) - 2.0*lgamma((double)(1+l) )
	+ lgamma((double)(2*l+n+2))
	)/(double)(l+n+1);
}


void CBrock::initialize(void)
{
  // Do nothing
}

void CBrock::get_dpotl(int lmax, int nmax, double r,
		       Eigen::MatrixXd& p, Eigen::MatrixXd& dp, int tid)
{
  double x  = r_to_xi(r);
  double dx = d_r_to_xi(r);

  double fac   = 0.5*sqrt(1.0 - x*x);
  double rfac  = sqrt(0.5*(1.0 - x));
  double drfac = -0.5/(1.0 - x*x);
  
  for (int l=0; l<=lmax; l++) {
    double dfac1 = 1.0 + x + 2.0*x*l;
    double dfac2 = 2.0*(l + 1);

    get_ultra(nmax-1, (double)l,     x, u[tid]);
    get_ultra(nmax-1, (double)(l+1), x, du[tid]);

    for (int n=0; n<nmax; n++) p(l, n) = rfac*u[tid][n];
    dp(l, 0) = dx*drfac*dfac1*rfac*u[tid][0];
    for (int n=1; n<nmax; n++) dp(l, n) = dx*rfac*(drfac*dfac1*u[tid][n] + 
						   dfac2*du[tid][n-1]);

    rfac *= fac;
  }

}

void CBrock::get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  double x    = r_to_xi(r);
  double fac  = 0.5*sqrt(1.0 - x*x);
  double rfac = sqrt(0.5*(1.0 - x));
  
  for (int l=0; l<=lmax; l++) {
    get_ultra(nmax-1, (double)l, x, u[tid]);
    for (int n=0; n<nmax; n++) p(l, n) = rfac*u[tid][n];
    rfac *= fac;
  }

}

void CBrock::get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  double x    = r_to_xi(r);
  double fac  = 0.5*sqrt(1.0 - x*x);
  double rfac = pow(0.5*(1.0 - x),2.5);

  for (int l=0; l<=lmax; l++) {
    get_ultra(nmax-1, (double)l, x, u[tid]);
    for (int n=0; n<nmax; n++) p(l, n) = krnl(l, n)*rfac*u[tid][n];
    rfac *= fac;
  }

}

void CBrock::get_potl_dens(int lmax, int nmax, double r, 
			   Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid)
{
  double x   = r_to_xi(r);
  double fac = 0.5*sqrt(1.0 - x*x);

  double rfac_p = sqrt(0.5*(1.0 - x));
  double rfac_d = pow(0.5*(1.0 - x),2.5);
  
  for (int l=0; l<=lmax; l++) {
    get_ultra(nmax-1, (double)l, x, u[tid]);
    for (int n=0; n<nmax; n++) {
      p(l, n) = rfac_p*u[tid][n];
      d(l, n) = krnl(l, n)*rfac_d*u[tid][n];
    }
    rfac_p *= fac;
    rfac_d *= fac;
  }
}


/*------------------------------------------------------------------------
 *                                                                       *
 *      Convert between reduced coordinate                               *
 *                                                                       *
 *               r^2-1                                                   *
 *          x = -------                                                  *
 *               r^2+1                                                   *
 *                                                                       *
 *      and its inverse:                                                 *
 *                                                                       *
 *              (1+x)^(1/2)                                              *
 *          r = ------------                                             *
 *              (1-x)^(1/2)                                              *
 *                                                                       *
 *-----------------------------------------------------------------------*/


#define BIG 1.0e30
double CBrock::xi_to_r(double xi)
{
  if (xi>=1.0) 
    return BIG;
  else
    return sqrt( (1.0+xi)/(1.0-xi) );
}

double CBrock::d_r_to_xi(double r)
{
  double fac;

  fac = r*r + 1.0;;
  return 4.0*r/(fac*fac);
}

double CBrock::r_to_xi(double r)
{
  return (r*r-1.0)/(r*r+1.0);
}
#undef BIG
