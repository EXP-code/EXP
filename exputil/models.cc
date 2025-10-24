#include "phase.H"
#include "models.H"

Eigen::Vector3d Isothermal::force(Eigen::Vector3d &u)
{
  double r, fr, x;
  static Eigen::Vector3d tmp;
  
  r = sqrt(u.dot(u));
  
  if (r <= TINY)
    {
      tmp.setZero();
      return tmp;
    }
  
  x = r/r0;
  
  fr =  -vcsq*r0*(x - atan(x))/r/r;
  
  tmp = fr*u/r;
  
  return tmp;
}

double Isothermal::potential(Eigen::Vector3d &u)
{
  double x, r;
  
  r = sqrt(u.dot(u));
  x = r/r0;
  
  return vcsq*(atan(x)/x + log(1.0+x*x)/2.0);
}




Eigen::Vector3d Miyamoto::force(Eigen::Vector3d &u)
{
  double rc, z, s, t, w;
  static Eigen::Vector3d tmp;
  
  tmp = u;
  
  rc = sqrt(SQR(u[1]) + SQR(u[2]));	/* cylindrical radius */
  z = u[3];
  
  t = sqrt(z*z + b*b);
  s = a + t;
  w = sqrt(s*s + rc*rc);
  w = w*w*w;
  
  tmp *= -GM/w;
  tmp[3] *= s/t;
  
  return tmp;
}

double Miyamoto::potential(Eigen::Vector3d &u)
{
  double rc, z, s, t, w;
  
  rc = sqrt(SQR(u[1]) + SQR(u[2]));	/* cylindrical radius */
  z = u[3];
  t = sqrt(z*z + b*b);
  s = a + t;
  w = sqrt(s*s + rc*rc);
  
  return -GM/w;
}


Eigen::Vector3d Hernquist::force(Eigen::Vector3d &u)
{
  double r, fr;
  static Eigen::Vector3d tmp;
  
  r = sqrt(u.dot(u));
  
  fr = -GM/SQR(r+r0);
  
  tmp = fr*u/r;
  return tmp;
}


double Hernquist::potential(Eigen::Vector3d &u)
{
  double r;
  r = sqrt(u.dot(u));
  
  return -GM/(r + r0);
}


double Logarithmic::potential(Eigen::Vector3d &u)
{
  return 0.5*vcsq*log(r0*r0 + u[1]*u[1] + SQR(u[2]/qy) 
		      + SQR(u[3]/qz));
}

Eigen::Vector3d Logarithmic::force(Eigen::Vector3d &u)
{
  double s;
  static Eigen::Vector3d tmp;
  tmp = u;
  
  s = r0*r0 + u[1]*u[1] + SQR(u[2]/qy) + SQR(u[3]/qz);
  
  tmp *= -vcsq/s;
  tmp[2] *= 1.0/qy/qy;
  tmp[3] *= 1.0/qz/qz;
  
  return tmp;
}



Eigen::Vector3d Plummer::force(Eigen::Vector3d &x)
{
  static Eigen::Vector3d f;
  double s;
  
  s = sqrt(1.0 + x.dot(x)/r2);
  
  f = -GM*x/s/s/s/r3;
  
  return f;
}

double Plummer::potential(Eigen::Vector3d &x)
{
  return -GM/r0/sqrt(1.0 + x.dot(x)/r2);
}
