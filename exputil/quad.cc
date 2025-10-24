
#include <cmath>
#include "numerical.H"
#include "models.H"


double Quadrupole_Bar::potential(Eigen::Vector3d &x)
{
  double s2;
  
  s2 = a2 + x.dot(x);
  
  return -0.5*GQ*(3.0*x[0]*x[0] - x.dot(x))/s2/s2/sqrt(s2);
}


Eigen::Vector3d Quadrupole_Bar::force(Eigen::Vector3d &x)
{
  static Eigen::Vector3d f;
  double r2, s2, s5, s7;
  
  r2 = x.dot(x);
  s2 = r2 + a2;
  s5 = s2*s2*sqrt(s2);
  s7 = s5*s2;
  
  f = 0.5*GQ*x/s7 * (3.0*r2 - 15.0*x[0]*x[0] - 2.0*a2);
  f[0] += 3.0*GQ*x[0]/s5;
  
  return f;
}
