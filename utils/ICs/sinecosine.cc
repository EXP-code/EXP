#include <cmath>
#include <Eigen/Dense>
  
void sinecosine_R(int mmax, double phi, Eigen::VectorXd& c, Eigen::VectorXd& s)
{

  c[0] = 1.0;
  s[0] = 0.0;

  c[1] = cos(phi);
  s[1] = sin(phi);

  for (int m=2; m<=mmax; m++) {
    c[m] = 2.0*c[1]*c[m-1] - c[m-2];
    s[m] = 2.0*c[1]*s[m-1] - s[m-2];
  }
}
