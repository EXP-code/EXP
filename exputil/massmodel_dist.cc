/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine returns a distribution function with desired anisotropy 
 *  according radius (radial, Merritt's Type I and circular, Merritt's Type
 *  II for r_max < ra)
 *
 *  Call sequence:
 *  -------------
 *
 *  oskipkov(ra,num);
 *
 *  void osipkov();
 *  double ra;
 *  int num;
 *
 *  Parameters:
 *  ----------
 *
 *  ra       anisotropy radius in units of R
 *  num      number of points in grid
 *
 *  Returns:
 *  -------
 *
 *  none
 *
 *  Notes:
 *  -----
 *
 *  Assumes that model does not have a density cusp, that is dRho/dr-->0 at
 *  r=0.  Third derivative boundary conditions (-1.0e30 input) are used for
 *  cubic splines.
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *  updated for Type II 2/24/89
 *  updated to deal with DIVERGE density types more accurately
 *
 ***************************************************************************/

#include <math.h>
#include <string>
#include <massmodel.H>
#include <interp.H>

#define OFFSET 1.0e-3
#define OFFTOL 1.0e-5

extern double gint_0(double a, double b, double (*f) (double), int NGauss);
extern double gint_2(double a, double b, double (*f) (double), int NGauss);

#define TSTEP 1.0e-8
#define NGauss 96

static int DIVERGE=0;

Eigen::VectorXd rhoQx;
Eigen::VectorXd rhoQy;
Eigen::VectorXd rhoQy2;
	     
void SphericalModelTable::setup_df(int NUM, double RA)
{
  double x,fac,fint(double p),d,dQ,Q,Qmin,Qmax;

  DIVERGE = diverge;

  df.ra2 = RA > 0.0 ? RA*RA : -RA*RA;
  if (df.ra2 < 0 && -RA < get_max_radius())
    bomb("Illegal value for osipkov radius");

/* Compute rho_Q(phi) */

  rhoQx.resize(num);
  rhoQy.resize(num);
  rhoQy2.resize(num);

  for (int i=0; i<num; i++) {
    x = density.x[i];
    rhoQx[i] = pot.y[i];
    if (diverge) {
      rhoQy[i] = (1.0 + x*x/df.ra2)*density.y[i]*pow(x,-diverge_rfac);
      rhoQy[i] = log(rhoQy[i] + TSTEP);
    }
    else
      rhoQy[i] = (1.0 + x*x/df.ra2)*density.y[i];
  }

  Spline(rhoQx, rhoQy, -1.0e30, -1.0e30, rhoQy2);
  

/* Tabulate the integral:

          Qmax
          /      1       d rho_Q
   I(Q) = | dp ------    -------
          /         1/2    dp
	  Q    [p-Q]

*/

  df.Q   .resize(NUM);
  df.fQ  .resize(NUM);
  df.ffQ .resize(NUM);
  df.fQ2 .resize(NUM);
  df.ffQ2.resize(NUM);
  df.num  = NUM;

  Qmax = pot.y[pot.num-1];
  Qmin = pot.y[0];
  dQ = (Qmax-Qmin)/(double)(df.num-1);
  
  df.Q[df.num-1] = Qmax;
  df.ffQ[df.num-1] = 0.0;
  fac = 1.0/(sqrt(8.0)*M_PI*M_PI);
  for (int i=df.num-2; i>=0; i--) {
    df.Q[i] = df.Q[i+1] - dQ;
    Q = df.Q[i];
    df.ffQ[i] = fac * gint_2(Q, Qmax, fint, NGauss);
  }

  df.off = df.Q[0] * ( 1.0 + OFFSET );
  for (int i=0; i<df.num; i++) {
    df.Q[i] = log(df.Q[i] - df.off);
    df.ffQ[i] = log(-df.ffQ[i] + OFFTOL);
  }

  Spline(df.Q, df.ffQ, -1.0e30,-1.0e30, df.ffQ2);


  /* Tabulate the df! */

  for (int i=df.num-1; i>=0; i--)
    Splint2(df.Q, df.ffQ, df.ffQ2, df.Q[i], d, df.fQ[i]);

  Spline(df.Q, df.fQ, -1.0e30, -1.0e30, df.fQ2);

  dist_defined = true;

}



double fint(double p)
{
  double d,d1;
  
  Splint2(rhoQx, rhoQy, rhoQy2, p, d, d1);
  if (DIVERGE) {
    d = exp(d) - TSTEP;
    d1 *= d + TSTEP;
  }
  return d1;
}



/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes dF(Q)/dQ for Osipkov-Merritt models
 *
 *
 *  Call sequence:
 *  -------------
 *  fq = distf(Q);          value of F(Q)
 *  dfq = dfdq(Q);          value of dF(Q)/dQ
 *
 *  double fq,dfq,Q;
 *
 *  Parameters:
 *  ----------
 *
 *  Q        Osipkov-Merritt Q = E + J^2/(2*r_a^2)
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *  Osipkov-Merritt model must be initialized by a call to osipkov();
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *
 ***************************************************************************/


double SphericalModelTable::distf(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, g, Q = log(E + 0.5*L*L/df.ra2 - df.off);

  if (Q >= df.Q[df.num-1])
    d = 0.0;
  else {
    Q = std::max<double>(df.Q[0], Q);
    Splint1(df.Q, df.fQ, df.fQ2, Q, d);
    Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
    d *= - exp(g - Q);
  }

  return d;
}


double SphericalModelTable::dfde(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, g, h, d1, Q = log(E + 0.5*L*L/df.ra2 - df.off);

  if (Q >= df.Q[df.num-1])
    d1 = 0.0;
  else {
    Q = std::max<double>(df.Q[0], Q);
    Splint2(df.Q, df.fQ, df.fQ2, Q, d, h);
    Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
    d1  = - exp(g - 2.0*Q) * ( h + d*(d - 1.0) );
  }

  return d1;
}


double SphericalModelTable::dfdl(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, g, h, d1, Q = log(E + 0.5*L*L/df.ra2 - df.off);

  if (Q >= df.Q[df.num-1])
    d1 = 0.0;
  else {
    Q = std::max<double>(df.Q[0], Q);
    Splint2(df.Q, df.fQ, df.fQ2, Q, d, h);
    Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
    d1  = - exp(g - 2.0*Q) * ( h + d*(d - 1.0) );
  }

  return d1*L/df.ra2;
}


double SphericalModelTable::d2fde2(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, g, h, k, d2, Q = E + 0.5*L*L/df.ra2;

  if (Q >= df.Q[df.num-1])
    d2 = 0.0;
  else {
    Q = std::max<double>(df.Q[0], Q);
    Splint3(df.Q, df.fQ, df.fQ2, Q, d, h, k);
    Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
    d2  = - exp(g - 3.0*Q) * (
			      (d - 2.0) * ( h + d*(d - 1.0) ) +
			      (k + 2.0*h*d - h)
			      );
  }

  return d2;
}

