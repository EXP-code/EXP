
#include <iostream>
#include <cmath>

#include "numerical.H"

#include "phase.H"
#include "staeckel.H"
#include "models.H"

/* 
   get the force at a 3D position *x.
*/


Eigen::Vector3d Perfect::force(Eigen::Vector3d &x)
{
  double r, fx, fr;
  static Eigen::Vector3d f;
  
  r = sqrt(SQR(x[2]) + SQR(x[3]));
  
  meridional_force(x[1], r, &fx, &fr);
  
  f[1] = fx;
  if (r >= 1.e-6)
    {
      f[2] = fr * x[2] / r;
      f[3] = fr * x[3] / r;
    }
  else
    {
      f[2] = 0.0;
      f[3] = 0.0;
    }
  
  return f;
}







/*
  return potential given 3D position 
*/



double Perfect::potential(Eigen::Vector3d &x)
{
  double r, lambda, mu, V;
  
  /* get the ellipsoidal coordinates */
  
  r = sqrt(x[1]*x[1] + x[2]*x[2]);
  to_ellipsoidal(x[0], r, &lambda, &mu);
  
  
  
  
  /* evaluate the potential */
  
  V = -1.0/(lambda-mu) * ((lambda + alpha)*Flambda(lambda) 
			  - (mu + alpha)*Fmu(mu));
  
  return V;
}




/*
  compute integrals of motion given Cartesian phase-space point
*/

void Perfect::Integrals(Phase &u, double *e, double *i2, double *jx)
{
  double r, vr;
  double vx, vy, x, y, z, vz;
  
  
  x = u.Position()[1];
  y = u.Position()[2];
  z = u.Position()[3];
  vx = u.Velocity()[1];
  vy = u.Velocity()[2];
  vz = u.Velocity()[3];
  r = sqrt(y*y + z*z);
  if (r >= 1.0e-12) vr = (y*vy + z*vz)/r;
  else vr = sqrt(vy*vy + vz*vz);
  
  
  *jx = y*vz - z*vy;
  
  *e = 0.5*(vx*vx + vy*vy + vz*vz) + potential(u.Position());
  
  *i2 = second_integral(x, r, vx, vr, *jx) 
    - alpha*(*e)  + SQR(*jx)/2.0  ;
  *i2 *= -1.0;
}




/* get the force in the x-R plane */



void Perfect::meridional_force(double x, double r, double *fx, double *fr)
{
  double lambda, mu, dldx, dldr, dmdr, dmdx, det;
  double Fl, Fm;
  double mu_alpha, lambda_alpha, dVmu, dVlambda;
  
  
  
  /* get the ellipsoidal coordinates */
  
  to_ellipsoidal(x, r, &lambda, &mu);
  
  
  
  
  // if determinant is zero, you've hit a marginal orbit. Sorry.
  
  det = lambda - mu;
  
  if (fabs(det) <= TINY) {
    std::cerr << "You have hit one of the focal points!" << std::endl;
    exit(-1);
  }
  
  
  // evaluate jacobian matrix of coordinate transformation
  
  dldx = 2.0*x*(beta + lambda)/det;
  dmdx = -2.0*x*(beta + mu)/det;
  dldr = 2.0*r*(alpha + lambda)/det;
  dmdr = -2.0*r*(alpha + mu)/det;
  
  
  
  
  
  /* a couple of handy auxiliary variables */
  
  mu_alpha = mu + alpha;
  lambda_alpha = lambda + alpha;
  
  
  Fl = Flambda(lambda);
  Fm = Fmu(mu);
  
  
  
  
  
  
  /* the derivatives of potential in ellipsoidal coords. */
  
  dVlambda = (lambda_alpha * Fl - mu_alpha*Fm)/det/det 
    - (Fl + lambda_alpha*dFlambda(lambda))/det;
  
  dVmu = -(lambda_alpha * Fl - mu_alpha * Fm) / det/det 
    + (Fm + mu_alpha*dFmu(mu))/det;
  
  
  
  
  /* now convert to cartesian components */
  
  *fx = -dVlambda * dldx - dVmu * dmdx;
  *fr = -dVlambda * dldr - dVmu * dmdr;
}	







/*
  get the ellipsoidal coordinates given the cartesian position
*/



void Perfect::to_ellipsoidal(double x, double r, double *lambda, double *mu)
{
  double x2, r2, tmp, q, b, c, delta;
  
  x2 = x*x;
  r2 = r*r;
  
  
  /* assemble coefficients of quadratic equation for tau */
  
  b = alpha + beta - x2 - r2;
  /* b is always negative */
  c = alpha*beta - x2*beta - r2*alpha;
  /* c is always positive */
  
  
  /* solve the quadratic */
  
  /* I don't check for positivity here, but have had no trouble */
  delta = sqrt(b*b - 4.0*c);
  q = -0.5*(b - delta);
  *lambda = q;
  *mu = c/q;
  
  
  
  
  /* sort the solutions */
  
  if (*mu > *lambda)
    {
      tmp = *mu;
      *mu = *lambda;
      *lambda = tmp;
    }
}



/* effective potential, with motion about the symmetry axis
   absorbed 
*/

double Perfect::Feff_lambda(double lambda, double jx)
{
  double I3, centrifugal;
  
  I3 = (beta - alpha) * jx*jx / 2.0;
  if (fabs(I3) >= TINY)
    {
      centrifugal = I3/(lambda + beta);
    }
  else
    {
      centrifugal = 0.0;
    }	
  
  
  return (lambda+alpha)*Flambda(lambda) + centrifugal;
}



double Perfect::Feff_mu(double mu, double jx)
{
  double I3, centrifugal;
  
  I3 = (beta - alpha) * jx*jx / 2.0;
  if (fabs(I3) >= TINY)
    {
      centrifugal = I3/(mu + beta);
    }
  else
    {
      centrifugal = 0.0;
    }
  
  return (mu+alpha)*Fmu(mu) + centrifugal;
}

/*
  Compute I2 from meridional position and angular momentum
*/


double Perfect::second_integral(double x, double r, double vx, 
				double vr, double jx)
{
  double lambda, mu, klam, kmu, Qlam, Qmu;
  double al, bl, am, bm, ab;
  
  
  /* get ellipsoidal coordinates */
  
  to_ellipsoidal(x, r, &lambda, &mu);
  
  
  
  /* define some auxiliary quantities */
  
  al = alpha + lambda;
  bl = beta + lambda;
  am = alpha + mu;
  bm = beta + mu;
  ab = alpha - beta;
  
  /* the rest has to be done in a rather obscure way, to avoid
     singularities */
  
  klam = r*x*vr*vx + 0.5/ab * (vx*vx * am*bl - vr*vr * al*bm);
  kmu =  r*x*vr*vx + 0.5/ab * (vx*vx * al*bm - vr*vr * am*bl);
  
  
  Qlam = klam - Feff_lambda(lambda, jx);
  Qmu = kmu - Feff_mu(mu, jx);
  
  return (lambda*Qmu - mu*Qlam)/(lambda - mu);
}
