/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Return Euler angle matrix (ref. Goldstein)
 *
 *  Call sequence:
 *  -------------
 *  include <Vector.h>
 *
 *  Matrix return_euler(double PHI, double THETA, double PSI, int BODY);
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
 *  MDW 03/29/93
 *
 ***************************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <Vector.h>

// #define TEST for test program

Matrix return_euler(double PHI, double THETA, double PSI, int BODY);


Matrix return_euler(double PHI, double THETA, double PSI, int BODY)
{
  double sph, cph, sth, cth, sps, cps;

  Matrix euler(1, 3, 1, 3);

  sph = sin(PHI);
  cph = cos(PHI);
  sth = sin(THETA);
  cth = cos(THETA);
  sps = sin(PSI);
  cps = cos(PSI);
  
  if (BODY) {

    euler[1][1] =  cps*cph - cth*sph*sps;
    euler[2][1] =  cps*sph + cth*cph*sps;
    euler[3][1] =  sps*sth;
    
    euler[1][2] = -sps*cph - cth*sph*cps;
    euler[2][2] = -sps*sph + cth*cph*cps;
    euler[3][2] =  cps*sth;
  
    euler[1][3] =  sth*sph;
    euler[2][3] = -sth*cph;
    euler[3][3] =  cth;

  }
  else {
    
    euler[1][1] =  cps*cph - cth*sph*sps;
    euler[1][2] =  cps*sph + cth*cph*sps;
    euler[1][3] =  sps*sth;
      
    euler[2][1] = -sps*cph - cth*sph*cps;
    euler[2][2] = -sps*sph + cth*cph*cps;
    euler[2][3] =  cps*sth;
      
    euler[3][1] =  sth*sph;
    euler[3][2] = -sth*cph;
    euler[3][3] =  cth;
    
  }

  return euler;
}

#ifdef TEST

main(int argc, char **argv)
{
  double phi, theta, psi;

  cout << "Phi, Theta, Psi: ";
  cin >> phi;
  cin >> theta;
  cin >> psi;

  Matrix euler0 = return_euler(phi, theta, psi, 0);
  Matrix euler1 = return_euler(phi, theta, psi, 1);

  cout << endl;
  euler0.print(cout);
  cout << endl;
  euler1.print(cout);

  Matrix eulert = euler0*euler1;
  cout << endl;
  eulert.print(cout);

}

#endif // TEST
