#include <cmath>

#include <Eigen/Eigen>

#include "numerical.H"
#include "models.H"


double Needle::potential(Eigen::Vector3d &x)
{
  double Sp, Sm, p;
  
  p = b*b + SQR(x[1]) + SQR(x[2]);
  Sp = sqrt(p + SQR(x[0]+a));
  Sm = sqrt(p + SQR(x[0]-a));
  
  if (a < TINY)
    {
      return -G*M/sqrt(p + SQR(x[1]));
    }
  
  return k*log((x[0]-a+Sm)/(x[0]+a+Sp));
}


Eigen::Vector3d Needle::force(Eigen::Vector3d &x)
{
  static Eigen::Vector3d f;
  double Sp, Sm, p, q, Sprod, Ssum;
  
  p = b*b + SQR(x[1]) + SQR(x[2]);
  Sp = sqrt(p + SQR(x[0]+a));
  Sm = sqrt(p + SQR(x[0]-a));
  Sprod = Sp*Sm;
  Ssum = Sp+Sm;
  q = kf/p/Sprod * (Ssum - 4.0*x[0]*x[0]/Ssum);
  
  
  f[1] = -4.0*kf*x[0]/Sprod/Ssum;
  f[2] = -q*x[1];
  f[3] = -q*x[2];
  
  return f;
}



double Miyamoto_Needle::potential(Eigen::Vector3d &x)
{
  double Tp, Tm, p, q;
  
  p = b + sqrt(c*c + x[2]*x[2]);
  q = SQR(x[1]) + p*p;
  Tp = sqrt(q + SQR(x[0]+a));
  Tm = sqrt(q + SQR(x[0]-a));
  
  if (a < TINY)
    {
      return -G*M/sqrt(q + SQR(x[1]));
    }
  
  return k*log((x[0]-a+Tm)/(x[0]+a+Tp));
}


Eigen::Vector3d Miyamoto_Needle::force(Eigen::Vector3d &x)
{
  static Eigen::Vector3d f;
  double Tp, Tm, p, q, w, Tprod, Tsum;
  
  p = b + sqrt(c*c + x[2]*x[2]);
  q = SQR(x[1]) + p*p;
  Tp = sqrt(q + SQR(x[0]+a));
  Tm = sqrt(q + SQR(x[0]-a));
  Tprod = Tp*Tm;
  Tsum = Tp+Tm;
  w = kf/q/Tprod * (Tsum - 4.0*x[0]*x[0]/Tsum);
  
  
  f[1] = -4.0*kf*x[0]/Tprod/Tsum;
  f[2] = -w*x[1];
  f[3] = -w*x[2]*p/(p-b);
  
  return f;
}


double Miyamoto_Needle::density(Eigen::Vector3d &x)
{
  double ap, am, B2;
  
  B2 = b*b+x[1]*x[1];
  ap = a + x[0];
  am = a - x[0];
  
  if (a>0.0) 
    {
      return M*b/4.0/Pi/a/B2 * (
				ap/sqrt(ap*ap + B2) + am/sqrt(am*am + B2)); 
    }
  else
    {
      return M*b/2.0/Pi/pow((B2+ap*ap), 1.5);
    }
}
