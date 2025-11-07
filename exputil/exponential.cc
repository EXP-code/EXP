//
// 2d exponetial disk
//

#include <cmath>
#include <cstdlib>

#include "massmodel.H"
#include "exponential.H"
#include "EXPmath.H"

double ExponentialDisk::get_pot(const double r)
{
  double y=0.5*r/a;

  double I0=EXPmath::cyl_bessel_i(0, y);
  double I1=EXPmath::cyl_bessel_i(1, y);
  double K0=EXPmath::cyl_bessel_k(0, y);
  double K1=EXPmath::cyl_bessel_k(1, y);

  return -M_PI*den0*r*(I0*K1 - I1*K0);
}
							

double ExponentialDisk::get_dpot(const double r)
{
  double y=0.5*r/a;

  double I0=EXPmath::cyl_bessel_i(0, y);
  double I1=EXPmath::cyl_bessel_i(1, y);
  double K0=EXPmath::cyl_bessel_k(0, y);
  double K1=EXPmath::cyl_bessel_k(1, y);

  return 2.0*M_PI*den0*y*(I0*K0 - I1*K1);
}

double ExponentialDisk::get_dpot2(const double r)
{
  double y=0.5*r/a;

  double I0=EXPmath::cyl_bessel_i(0, y);
  double I1=EXPmath::cyl_bessel_i(1, y);
  double K0=EXPmath::cyl_bessel_k(0, y);
  double K1=EXPmath::cyl_bessel_k(1, y);

  return M_PI*den0/a*(I0*K0 + I1*K1 - 2.0*(I1*K0 - I0*K1));
}

void ExponentialDisk::get_pot_dpot(const double r, double &ur, double &dur)
{
  double y=0.5*r/a;
  
  double I0=EXPmath::cyl_bessel_i(0, y);
  double I1=EXPmath::cyl_bessel_i(1, y);
  double K0=EXPmath::cyl_bessel_k(0, y);
  double K1=EXPmath::cyl_bessel_k(1, y);

  ur = -M_PI*den0*r*(I0*K1 - I1*K0);
  dur = 2.0*M_PI*den0*y*(I0*K0 - I1*K1);
}
  
