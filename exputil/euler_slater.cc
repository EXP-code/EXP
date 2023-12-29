/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Return Euler angle matrix (ref. Slater)
 *
 *  Call sequence:
 *  -------------
 *
 *  Eigen::Matrix3d return_euler(double PHI, double THETA, double PSI, int BODY);
 *
 *  Euler angles: PHI, THETA, PSI
 *
 *  BODY = 0:     rotate axes, keep vector (or "body") fixed in space
 *  BODY = 1:     rotate body
 *
 *
 *  Parameters:
 *  ----------
 *
 *  TEST -- compile with main() test routine
 *
 *  Returns:
 *  -------
 *
 *  Euler transformation
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 09/6/96
 *
 ***************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <Eigen/Eigen>

// #define TEST for test program

Eigen::Matrix3d return_euler_slater(double PHI, double THETA, double PSI, int BODY)
{
  double sph, cph, sth, cth, sps, cps;

  Eigen::Matrix3d euler;

  sph = sin(PHI);
  cph = cos(PHI);
  sth = sin(THETA);
  cth = cos(THETA);
  sps = sin(PSI);
  cps = cos(PSI);
  
  euler(0, 0) = -sps*sph + cth*cph*cps;
  euler(0, 1) =  sps*cph + cth*sph*cps;
  euler(0, 2) =  cps*sth;
      
  euler(1, 0) = -cps*sph - cth*cph*sps;
  euler(1, 1) =  cps*cph - cth*sph*sps;
  euler(1, 2) = -sps*sth;
      
  euler(2, 0) = -sth*cph;
  euler(2, 1) = -sth*sph;
  euler(2, 2) =  cth;

  if (BODY)
    return euler.transpose();
  else
    return euler;
}

#ifdef TEST

int main(int argc, char **argv)
{
  const double onedeg = M_PI/180.0;
  double phi, theta, psi;

  std::cout << "Phi, Theta, Psi (in degrees): ";
  std::cin >> phi;
  std::cin >> theta;
  std::cin >> psi;

  phi   *= onedeg;
  theta *= onedeg;
  psi   *= onedeg;

  auto euler0 = return_euler_slater(phi, theta, psi, 0);
  auto euler1 = return_euler_slater(phi, theta, psi, 1);

  std::cout << std::endl << euler0 << std::endl;
  std::cout << std::endl << euler1 << std::endl;

  auto eulert = euler0*euler1;
  std::cout << std::endl << eulert << std::endl;

  return 0;
}

#endif // TEST
