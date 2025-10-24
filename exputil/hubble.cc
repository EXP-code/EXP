#include <cmath>
#include "numerical.H"
#include "models.H"


/* modified Hubble model */


double Modified_Hubble::potential(Eigen::Vector3d &x)
{
  double r, s;
  
  r = sqrt(x.dot(x));
  s = r/rs;
  if (r < TINY)
    {
      return -Gravity*Ms/rs;
    }
  
  return -Gravity*Ms/r * log(s + sqrt(s*s + 1.0));
}

Eigen::Vector3d Modified_Hubble::force(Eigen::Vector3d &x)
{
  static Eigen::Vector3d f;
  double r, s, fr;
  
  r = sqrt(x.dot(x));
  s = r/rs;
  if (r < TINY)
    {
      f.setZero();
      return f;
    }
  
  fr = Gravity*Ms/r/r * (r / sqrt(r*r + rs*rs) - 
			 log(s + sqrt(1.0 + s*s)));
  f = x*fr/r;
  return f;
}
