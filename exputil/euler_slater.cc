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

#include <cstdlib>
#include <cmath>
#include <Eigen/Eigen>

// #define TEST for test program

Eigen::MatrixXd return_euler_slater(double PHI, double THETA, double PSI, int BODY);


Eigen::MatrixXd return_euler_slater(double PHI, double THETA, double PSI, int BODY)
{
  double sph, cph, sth, cth, sps, cps;

  Eigen::MatrixXd euler(3, 3);

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

void main(int argc, char **argv)
{
  const double onedeg = M_PI/180.0;
  double phi, theta, psi;

  cout << "Phi, Theta, Psi (in degrees): ";
  cin >> phi;
  cin >> theta;
  cin >> psi;

  phi   *= onedeg;
  theta *= onedeg;
  psi   *= onedeg;

  Matrix euler0 = return_euler_slater(phi, theta, psi, 0);
  Matrix euler1 = return_euler_slater(phi, theta, psi, 1);

  cout << endl;
  euler0.print(cout);
  cout << endl;
  euler1.print(cout);

  Matrix eulert = euler0*euler1;
  cout << endl;
  eulert.print(cout);

}

#endif // TEST
