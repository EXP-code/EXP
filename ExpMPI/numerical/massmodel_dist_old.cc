// This may look like C code, but it is really -*- C++ -*-

const char rcsid_massmodel_dist[] = "$Id$";

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
#include <String.h>
#include <Vector.h>
#include <massmodel.h>
#include <iomanip.h>

void Spline(Vector &x, Vector &y, double yp1, double ypn, Vector &y2);
void Splint1(Vector &xa, Vector &ya, Vector &y2a, double x, double &y, int even=0);
void Splint2(Vector &xa, Vector &ya, Vector &y2a, double x, double &y, 
	     double &dy, int even=0);
void Splint3(Vector &xa, Vector &ya, Vector &y2a, double x, double &y, 
	     double &dy, double &ddy, int even=0);

extern "C" {
  double gint_0(double a, double b, double (*f) (double), int NGauss);
  double gint_2(double a, double b, double (*f) (double), int NGauss);
};

#define TSTEP 1.0e-8
#define NGauss 96

static int DIVERGE=0;

Vector rhoQx;
Vector rhoQy;
Vector rhoQy2;
	     
void SphericalModelTable::setup_df(int NUM, double RA)
{
  double x,ur,fac,fint(double p),d,d1,d2,dQ,Q,Qmin,Qmax;
  int i;

  DIVERGE = diverge;

  df.ra2 = RA > 0.0 ? RA*RA : -RA*RA;
  if (df.ra2 < 0 && -RA < get_max_radius())
    bomb("Illegal value for osipkov radius");

/* Compute rho_Q(phi) */

  rhoQx.setsize(1, num);
  rhoQy.setsize(1, num);
  rhoQy2.setsize(1, num);

  for (i=1; i<=num; i++) {
    x = density.x[i];
    rhoQx[i] = pot.y[i];
    if (diverge) {
      rhoQy[i] = (1.0 + x*x/df.ra2)*density.y[i]*pow(x,-diverge_rfac);
      rhoQy[i] = log(rhoQy[i] + TSTEP);
    }
    else
      rhoQy[i] = (1.0 + x*x/df.ra2)*density.y[i];
  }

  Spline(rhoQx, rhoQy, -1.0e30,-1.0e30, rhoQy2);
  

/* Tabulate the integral:

          Qmax
          /      1       d rho_Q
   I(Q) = | dp ------    -------
          /         1/2    dp
	  Q    [p-Q]

*/

  df.Q    = Vector(1, NUM);
  df.fQ   = Vector(1, NUM);
  df.fQ2  = Vector(1, NUM);
  Vector ffQ(1, NUM);
  Vector ffQ2(1, NUM);
  df.num  = NUM;

  Qmax = pot.y[pot.num];
  Qmin = pot.y[1];
  dQ = (Qmax-Qmin)/(double)(df.num-1);
  
  df.Q[df.num] = Qmax;
  ffQ[df.num] = 0.0;
  fac = 1.0/(sqrt(8.0)*M_PI*M_PI);
  for (i=df.num-1; i>0; i--) {
    df.Q[i] = df.Q[i+1] - dQ;
    Q = df.Q[i];
    ffQ[i] = fac * gint_2(Q, Qmax, fint, NGauss);
  }
  Spline(df.Q, ffQ, -1.0e30,-1.0e30, ffQ2);


/* Tabulate the df! */

  for (i=df.num; i>0; i--) {
    Splint2(df.Q, ffQ, ffQ2, df.Q[i], d, df.fQ[i]);
    if (df.fQ[i]<0.0) {
      cerr << "SphericalModelTable : distribution has negative density!\n";
#ifdef DEBUG
      cerr << "SphericalModelTable : Q=" << df.Q[i] <<"  F=" << df.fQ[i] << endl;
#else
      _exit(-1);
#endif
    }
  }

  Spline(df.Q, df.fQ, -1.0e30, -1.0e30, df.fQ2);

  dist_defined = TRUE;

}



double fint(double p)
{
  double d,d1,d2;
  
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

  double d, Q = E + 0.5*L*L/df.ra2;

  if (Q > df.Q[df.num])
    d = 0.0;
  else
    Splint1(df.Q, df.fQ, df.fQ2, Q, d);

  return d;
}


double SphericalModelTable::dfde(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, d1, Q = E + 0.5*L*L/df.ra2;

  if (Q > df.Q[df.num])
    d1 = 0.0;
  else
    Splint2(df.Q, df.fQ, df.fQ2, Q, d, d1);

  return d1;
}


double SphericalModelTable::dfdl(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, d1, Q = E + 0.5*L*L/df.ra2;

  if (Q > df.Q[df.num])
    d1 = 0.0;
  else
    Splint2(df.Q, df.fQ, df.fQ2, Q, d, d1);

  return d1*L/df.ra2;
}


double SphericalModelTable::d2fde2(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, d1, d2, Q = E + 0.5*L*L/df.ra2;

  if (Q > df.Q[df.num])
    d1 = 0.0;
  else
    Splint3(df.Q, df.fQ, df.fQ2, Q, d, d1, d2);

  return d2;
}


