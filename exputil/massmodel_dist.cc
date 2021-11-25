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
#define OFFTOL 1.2

extern double gint_0(double a, double b, double (*f) (double), int NGauss);
extern double gint_2(double a, double b, double (*f) (double), int NGauss);

#define TSTEP 1.0e-8
#define NGauss 96

static int DIVERGE=0;

Eigen::VectorXd rhoQx;
Eigen::VectorXd rhoQy;
Eigen::VectorXd rhoQy2;
	     
void SphericalModelTable::debug_fdist()
{
  static int cur_count = 0;
  std::ostringstream sout;
  sout << "test_fdist." << cur_count++;
  std::ostringstream out(sout.str());
  if (out) {
    if (chebyN) {
      for (int n=0; n<dfc.num; n++)
	out << std::setw(15) << dfc.Q[n]
	    << std::setw(15) << dfc.fQ[n]
	    << std::setw(15) << dfc.ffQ[n]
	    << std::setw(15) << dfc.FF.eval(dfc.Q[n])
	    << std::setw(15) << dfc.FF.deriv(dfc.Q[n])
	    << std::endl;

    } else {
      for (int n=0; n<df.num; n++)
	out << std::setw(15) << df.Q[n]
	    << std::setw(15) << df.fQ[n]
	    << std::setw(15) << df.ffQ[n]
	    << std::setw(15) << df.fQ2[n]
	    << std::setw(15) << df.ffQ2[n]
	    << std::endl;
    }
  }
}


void SphericalModelTable::setup_df(int NUM, double RA)
{
  double x,fac,fint(double p),d,dQ,Q,Qmin,Qmax;

  DIVERGE = diverge;

  double ra2 = RA > 0.0 ? RA*RA : -RA*RA;
  if (ra2 < 0 && -RA < get_max_radius())
    bomb("Illegal value for osipkov radius");

  // Compute rho_Q(phi)
  //
  rhoQx.resize(num);
  rhoQy.resize(num);
  rhoQy2.resize(num);

  for (int i=0; i<num; i++) {
    x = density.x[i];
    rhoQx[i] = pot.y[i];
    if (diverge) {
      rhoQy[i] = (1.0 + x*x/ra2)*density.y[i]*pow(x,-diverge_rfac);
      rhoQy[i] = log(rhoQy[i] + TSTEP);
    }
    else
      rhoQy[i] = (1.0 + x*x/ra2)*density.y[i];
  }

  Spline(rhoQx, rhoQy, -1.0e30, -1.0e30, rhoQy2);
  

/* Tabulate the integral:

          Qmax
          /      1       d rho_Q
   I(Q) = | dp ------    -------
          /         1/2    dp
	  Q    [p-Q]

*/

  if (chebyN) {

    dfc.Q  .resize(NUM);
    dfc.fQ .resize(NUM);
    dfc.ffQ.resize(NUM);
    dfc.num = NUM;
    dfc.ra2 = ra2;

    Qmax = get_pot(pot.x[pot.num-1]);
    Qmin = get_pot(pot.x[0]);
    dQ = (Qmax - Qmin)/(double)(dfc.num-1);
  
    foffset = -std::numeric_limits<double>::max();
    dfc.Q[dfc.num-1] = Qmax;
    dfc.ffQ[dfc.num-1] = 0.0;
    fac = 1.0/(sqrt(8.0)*M_PI*M_PI);
    for (int i=dfc.num-2; i>=0; i--) {
      dfc.Q[i] = dfc.Q[i+1] - dQ;
      Q = dfc.Q[i];
      dfc.ffQ[i] = fac * gint_2(Q, Qmax, fint, NGauss);
      foffset = std::max<double>(foffset, dfc.ffQ[i]);
    }
    
    // 'foffset' is now the largest value of the integral. Shift a bit
    // larger.
    //
    if      (foffset >  0.0) foffset *= OFFTOL;
    else if (foffset == 0.0) foffset = 1.0e-5;
    else                     foffset /= OFFTOL;

    dfc.off = dfc.Q[0] * ( 1.0 + OFFSET );
    for (int i=0; i<dfc.num; i++) {
      dfc.Q[i] = log(dfc.Q[i] - dfc.off);
      dfc.ffQ[i] = log(-dfc.ffQ[i] + foffset);
    }
    
    dfc.GG = Cheby1d(dfc.Q, dfc.ffQ, chebyN);
    
    for (int i=0; i<dfc.num; i++) {
      dfc.fQ[i] = dfc.GG.deriv(dfc.Q[i]);
    }

    dfc.FF = Cheby1d(dfc.Q, dfc.fQ, chebyN);

  } else {

    df.Q   .resize(NUM);
    df.fQ  .resize(NUM);
    df.ffQ .resize(NUM);
    df.fQ2 .resize(NUM);
    df.ffQ2.resize(NUM);
    df.num = NUM;
    df.ra2 = ra2;

    Qmax = pot.y[pot.num-1];
    Qmin = pot.y[0];
    dQ = (Qmax-Qmin)/(double)(df.num-1);
  
    df.Q[df.num-1] = Qmax;
    df.ffQ[df.num-1] = 0.0;
    fac = 1.0/(sqrt(8.0)*M_PI*M_PI);
    foffset = -std::numeric_limits<double>::max();
    for (int i=df.num-2; i>=0; i--) {
      df.Q[i] = df.Q[i+1] - dQ;
      Q = df.Q[i];
      df.ffQ[i] = fac * gint_2(Q, Qmax, fint, NGauss);
      foffset = std::max<double>(foffset, df.ffQ[i]);
    }
    
    // 'foffset' is now the largest value of the integral. Shift a bit
    // larger.
    //
    if      (foffset > 0.0)  foffset *= OFFTOL;
    else if (foffset == 0.0) foffset = 1.0e-5;
    else                     foffset /= OFFTOL;

    df.off = df.Q[0] * ( 1.0 + OFFSET );
    for (int i=0; i<df.num; i++) {
      df.Q[i] = log(df.Q[i] - df.off);
      df.ffQ[i] = log(-df.ffQ[i] + foffset);
    }
    
    Spline(df.Q, df.ffQ, -1.0e30,-1.0e30, df.ffQ2);

    
    // Tabulate the df!
    //
    if (linear) {
      for (int i=df.num-1; i>=0; i--)
	df.fQ[i] = drv2(df.Q[i], df.Q, df.ffQ);
    } else {
      for (int i=df.num-1; i>=0; i--)
	Splint2(df.Q, df.ffQ, df.ffQ2, df.Q[i], d, df.fQ[i]);
    
      Spline(df.Q, df.fQ, -1.0e30, -1.0e30, df.fQ2);
    }
  }

  dist_defined = true;

  debug_fdist();
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

  double d, g;

  if (chebyN) {

    double Q = E + 0.5*L*L/dfc.ra2 - dfc.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(dfc.Q[0], log(Q));

    if (Q > dfc.Q[dfc.num-1])
      d = 0.0;
    else {
      g = dfc.GG.eval(Q);
      d = dfc.FF.eval(Q);
      d *= - exp(g - Q);
    } 

  } else {
    
    double Q = E + 0.5*L*L/df.ra2 - df.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(df.Q[0], log(Q));

    if (Q > df.Q[df.num-1])
      d = 0.0;
    else {
      if (linear) {
	d = odd2(Q, df.Q, df.fQ);
	g = odd2(Q, df.Q, df.ffQ);
      } else {
	Splint1(df.Q, df.fQ, df.fQ2, Q, d);
	Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
      }
      d *= - exp(g - Q);
    }
  }

  return d;
}


double SphericalModelTable::dfde(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, g, h, d1;

  if (chebyN) {
    
    double Q = E + 0.5*L*L/dfc.ra2 - dfc.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(dfc.Q[0], log(Q));

    if (Q > dfc.Q[dfc.num])
      d1 = 0.0;
    else {
      d = dfc.FF.eval(Q);
      h = dfc.FF.deriv(Q);
      g = dfc.GG.eval(Q);
      d1 = - exp(g - 2.0*Q) * ( h + d*(d - 1.0) );
    }

  } else {

    double Q = E + 0.5*L*L/df.ra2 - df.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(df.Q[0], log(Q));

    if (Q > df.Q[df.num-1])
      d1 = 0.0;
    else {

      double Q = log(E + 0.5*L*L/df.ra2 - df.off);

      Q = max<double>(df.Q[0], Q);

      if (linear) {
	double dq = 1.0e-3*(df.Q[1]-df.Q[0]);
	d = odd2(Q, df.Q, df.fQ);
	h = (
	     odd2(Q+dq, df.Q, df.fQ) - 
	     odd2(Q-dq, df.Q, df.fQ)
	     ) / (2.0*dq);
	g = odd2(Q, df.Q, df.ffQ);
      } else {
	Splint2(df.Q, df.fQ, df.fQ2, Q, d, h);
	Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
      }
      d1  = - exp(g - 2.0*Q) * ( h + d*(d - 1.0) );
    }
  }
  
  return d1;
}


double SphericalModelTable::dfdl(double E, double L)
{

  if (!dist_defined) bomb("distribution function not defined");

  double d, g, h, d1;

  if (chebyN) {
    
    double Q = E + 0.5*L*L/dfc.ra2 - dfc.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(dfc.Q[0], log(Q));

    if (Q > dfc.Q[dfc.num])
      d1 = 0.0;
    else {
      d = dfc.FF.eval(Q);
      h = dfc.FF.deriv(Q);
      g = dfc.GG.eval(Q);
      d1 = - exp(g - 2.0*Q) * ( h + d*(d - 1.0) );
    }

  } else {

    double Q = E + 0.5*L*L/df.ra2 - df.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(df.Q[0], Q);

    if (Q > df.Q[df.num-1])
      d1 = 0.0;
    else {
      if (linear) {
	double dq = 1.0e-3*(df.Q[2]-df.Q[1]);
	d = odd2(Q, df.Q, df.fQ);
	h = (
	     odd2(Q+dq, df.Q, df.fQ) - 
	     odd2(Q-dq, df.Q, df.fQ)
	     ) / (2.0*dq);
	g = odd2(Q, df.Q, df.ffQ);
      } else {
	Splint2(df.Q, df.fQ, df.fQ2, Q, d, h);
	Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
      }
      d1  = - exp(g - 2.0*Q) * ( h + d*(d - 1.0) );
    }
  }

  return d1*L/df.ra2;
}


double SphericalModelTable::d2fde2(double E, double L)
{
  if (!dist_defined) bomb("distribution function not defined");

  double d, g, h, k, d2;

  if (chebyN) {

    double Q = E + 0.5*L*L/dfc.ra2 - dfc.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(dfc.Q[0], log(Q));
    
    if (Q > dfc.Q[dfc.num])
      d2 = 0.0;
    else {
      d = dfc.FF.eval(Q);
      h = dfc.FF.deriv(Q);
      k = dfc.FF.deriv2(Q);
      g = dfc.GG.eval(Q);

      d2  = - exp(g - 3.0*Q) * (
				(d - 2.0) * ( h + d*(d - 1.0) ) +
				(k + 2.0*h*d - h)
				);

    } 
    
  } else {

    double Q = E + 0.5*L*L/df.ra2 - df.off;
    if (Q<=0.0) return 0.0;

    Q = max<double>(df.Q[0], log(Q));

    if (Q > df.Q[df.num-1])
      d2 = 0.0;
    else {
      if (linear) {
	double dq = 1.0e-3*(df.Q[2]-df.Q[1]);
	d = odd2(Q, df.Q, df.fQ);
	h = (
	   odd2(Q+dq, df.Q, df.fQ) - 
	   odd2(Q-dq, df.Q, df.fQ)
	   ) / (2.0*dq);
	k = (
	     odd2(Q+2.0*dq, df.Q, df.fQ) -
	     2.0*d +
	     odd2(Q-2.0*dq, df.Q, df.fQ)
	     ) / (4.0*dq*dq);
	g = odd2(Q, df.Q, df.ffQ);
      } else {
	Splint3(df.Q, df.fQ, df.fQ2, Q, d, h, k);
	Splint1(df.Q, df.ffQ, df.ffQ2, Q, g);
	d2  = - exp(g - 3.0*Q) * (
				  (d - 2.0) * ( h + d*(d - 1.0) ) +
				  (k + 2.0*h*d - h)
				  );
      }
    }
  }

  return d2;
}

