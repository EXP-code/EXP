#include <iostream>
#include <cmath>
#include "numerical.H"
#include "rotcurv.H"
#include "phase.H"

/*
  ROTATION CURVE ANALYSIS FUNCTIONS
  
  double v_circ(r)
  double omega_circ(r)
  double epicyclic_kappa(r)
  double Oort_A(r)
  double Oort_B(r)
  
  These functions work with the axisymmetric component of the radial
  force. The force is azimuthally averaged by numerical integration,
  whether it needs it or not.
  
  A force function of the form
  Force3D(x, f)
  is assumed.
  
  Any derivatives are done numerically.
  
  
  
  
  Written by:
  KL 8/18/91
  ANSI Version 9/30/91
  3D Version 10/6/91
*/


extern Frame_Rotation Frame;

//	Eigen::Vector3d Force(double, Eigen::Vector3d &, Eigen::Vector3d &);







/* a variable with file scope */

static double global_rdum;


double v_circ(double r)
{
  double axisymmetric_force(double), vc, tmp, tmp2, tmpsq;
  
  tmp = Frame.omega;
  tmpsq = Frame.omegasq;
  tmp2 = Frame.omega2;
  Frame.omega = 0.0;
  Frame.omega2 = 0.0;
  Frame.omegasq = 0.0;
  
  vc = sqrt(r * fabs(axisymmetric_force(r)));
  
  Frame.omega = tmp;
  Frame.omegasq = tmpsq;
  Frame.omega2 = tmp2;
  
  return vc;
}

double omega_circ(double r)
{
  double h=1.e-4, axisymmetric_force(double);
  double tmp, tmpsq, tmp2, w;
  
  if (r > TINY) return v_circ(r)/r;
  
  /* special case at r=0: use F=-omega^2 r.   */
  
  tmp = Frame.omega;
  tmpsq = Frame.omegasq;
  tmp2 = Frame.omega2;
  Frame.omega = 0.0;
  Frame.omega2 = 0.0;
  Frame.omegasq = 0.0;
  
  w =  sqrt(fabs(axisymmetric_force(h)/h));
  
  Frame.omega = tmp;
  Frame.omegasq = tmpsq;
  Frame.omega2 = tmp2;
  
  return w;
}

double Oort_A(double r)
{
  double h=1.e-4, dwdr;
  
  if (r > h) dwdr = (omega_circ(r + h) - omega_circ(r - h)) / 2.0/h;
  else dwdr = (omega_circ(r + h) - omega_circ(r))/h;
  
  return -0.5*r*dwdr;
}


double Oort_B(double r)
{
  
  return Oort_A(r) - omega_circ(r);
}


double epicyclic_kappa(double r)
{
  return sqrt(fabs(4.0 * Oort_B(r) * omega_circ(r)));
}

double vertical_kappa(double r)
{
  double vert_fdum(double);
  
  global_rdum = r;
  
  /* average over azimuth */
  
  return qadapt(0.0, Pi/2.0, vert_fdum, 1.e-5)*2.0/Pi;
}	








/* get the axisymmetrix component of the radial force */


double axisymmetric_force(double r)
{
  double fdum(double);
  
  global_rdum = r;
  
  /* average over azimuth */
  
  return (fdum(0.0) + fdum(Pi/2.0))/2.0;
  
  /*
    return qadapt(0.0, Pi/2.0, fdum, 1.e-5)*2.0/Pi;
  */
}	




/* 
   The force as a function of azimuthal angle. The radial dependence
   is passed through the global variable global_rdum.
*/

double fdum(double theta)
{
  Eigen::Vector3d x, v, f;
  v.setZero();
  
  if (global_rdum <= TINY) return 0.0;
  
  x[0] = global_rdum * cos(theta);
  x[1] = global_rdum * sin(theta);
  x[2] = 0.0;
  
  f = (*Phase::Force)(0.0, x, v);
  
  return x.dot(f)/global_rdum;	
}


double vert_fdum(double theta)
{
  double eps= 1.e-2;
  Eigen::Vector3d x, v, f;
  v.setZero();
  
  
  x[0] = global_rdum * cos(theta);
  x[1] = global_rdum * sin(theta);
  x[2] = eps;
  
  // remove centrifugal force 
  
  f = (*Phase::Force)(0.0, x, v) - Frame.omega*x;
  return sqrt(fabs(f[3]/eps));
}
