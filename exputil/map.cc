#include <cmath>
#include "numerical.H"
#include "phase.H"

/* 
   Run an orbit from one point on the surface of section 
   to another. 
*/

Phase Phase::x_return_map(void)
{
  static Phase tmp, u0, u1;
  double dt, tcross, xcross, vcross, ycross, vycross, zcross, vzcross;
  double vystart;
  int i, first_step, no_intersection, opposite_direction;
  
  u1 = *this;
  dt = get_Xscale()/get_Vscale();
  vystart = u1.v[2];
  
  
  /* run the orbit until it again crosses the surface of section */
  
  i=0;
  do {
    u0 = u1;
    u1 = u0.advance(dt);
    first_step = i<1;
    no_intersection = u0.x[2]*u1.x[2] >= 0.0;
    opposite_direction = u1.v[2]*vystart < 0.0;
    i++;
  } while (first_step || no_intersection);
  
  
  
  
  
  /* do inverse cubic interpolation to get precise corssing time */
  
  tcross = inverse_interp_zero(u0.x[2], u1.x[2], u0.v[2], u1.v[2], 
			       u0.t, u1.t);
  
  
  
  /* now get the position and velocity at the crossing */
  
  forward_interp(u0.x[1], u1.x[1], u0.v[1], u1.v[1], u0.t, u1.t,
		 tcross, xcross, vcross);
  
  forward_interp(u0.x[2], u1.x[2], u0.v[2], u1.v[2], u0.t, u1.t, 
		 tcross, ycross, vycross);
  
  forward_interp(u0.x[3], u1.x[3], u0.v[3], u1.v[3], u0.t, u1.t, 
		 tcross, zcross, vzcross);
  
  tmp.x.setZero();
  tmp.v.setZero();
  tmp.t = 0.0;
  
  tmp.x[1] = xcross;
  tmp.x[2] = ycross;
  tmp.x[3] = zcross;
  tmp.v[1] = vcross;
  tmp.v[2] = vycross;
  tmp.v[3] = vzcross;
  
  return tmp;
}



/* various cubic interpolation routines */

static double apoly[4];

void get_cubic_coeffs(double f0, double f1, double df0, 
		      double df1, double t0, double t1)
{
  double p, q, dt;
  
  dt = t1 - t0;
  apoly[0] = f0;
  apoly[1] = df0*dt;
  
  p = f1 - f0 - apoly[1];
  q = (df1- df0)*dt;
  apoly[2] = 3.0*p - q;
  apoly[3] = p - apoly[2];
}

void cubic_interp(double u, double& f, double& df)
{
  double dpoly[3];
  
  dpoly[0]=apoly[1]; dpoly[1]=2.0*apoly[2]; dpoly[2]=3.0*apoly[3];
  f = apoly[0] + u*(apoly[1] + u*(apoly[2] + u*apoly[3]));
  df = dpoly[0] + u*(dpoly[1] + u*dpoly[2]);
}



double eval_interp(double u)
{
  return poly(u, 3, apoly);
}

double inverse_interp_zero(double f0, double f1, double df0,
			   double df1, double t0, double t1)
{
  get_cubic_coeffs(f0, f1, df0, df1, t0, t1);

  double u = rtsafe(cubic_interp, 0.0, 1.0, 1.e-5);
  double t = t0 + (t1-t0)*u;
  
  return t;
}

double inverse_interp_min(double f0, double f1, double df0, 
			  double df1, double t0, double t1)
{
  double u, t, fmin, eval_interp(double);
  
  get_cubic_coeffs(f0, f1, df0, df1, t0, t1);
  u = brent(0.0, 0.5, 1.0, eval_interp, 1.0e-5, fmin);
  t = t0 + (t1-t0)*u;
  
  return t;
}

void forward_interp(double f0, double f1, double df0, 
		    double df1, double t0, double t1, double t,
		    double& f, double& df)
{
  double u = (t-t0)/(t1-t0);
  get_cubic_coeffs(f0, f1, df0, df1, t0, t1);
  cubic_interp(u, f, df);
  df = df / (t1-t0);
}
