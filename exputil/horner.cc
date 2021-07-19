/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine returns the results of a Horner's scheme for a Poly
 *  Ref: Henrici, Applied and Computational Complex Analysis, p. 433f
 *
 *  Call sequence:
 *  -------------
 *
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 01/07/93
 *
 ***************************************************************************/

#include <cmath>
#include <cstdlib>
#include <complex>
#include <Eigen/Eigen>

#include "cpoly.H"

Eigen::VectorXd get_horner(double z, Poly &p)
{
  int order = p.getorder();

  Eigen::MatrixXd b(order+2, order+1);
  for (int i=0; i<=order; i++) b(0, i) = p[order-i];
  
  for (int i=0; i<=order; i++) {
    b(i+1, 0) = b(0, 0);
    for (int j=1; j<=order-i; j++)
      b(i+1, j) = b(i, j) + z*b(i+1, j-1);
  }

  Eigen::VectorXd ans(order+1);
  for (int i=0; i<=order; i++) ans[i] = b(i+1, order-i);
  
  return ans;
}


Eigen::VectorXcd Cget_horner(std::complex<double> z, CPoly &p)
{
  int order = p.getorder();
  Eigen::VectorXcd ans;

  if (order==0) {
    ans.resize(2);
    ans[0] = p[0];
    ans[1] = 0.0;
    return ans;
  }

  Eigen::MatrixXcd b(order+2, order+1);
  for (int i=0; i<=order; i++) b(0, i) = p[order-i];
  
  for (int i=0; i<=order; i++) {
    b(i+1, 0) = b(0, 0);
    for (int j=1; j<=order-i; j++)
      b(i+1, j) = b(i, j) + z*b(i+1, j-1);
  }

  ans.resize(order+1);
  for (int i=0; i<=order; i++) ans[i] = b(i+1, order-i);
  
  return ans;
}

