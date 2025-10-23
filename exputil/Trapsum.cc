#include <cstdlib>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <cmath>

#include "interp.H"


double Trapsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
{
  int sz = x.size();

  double h, p=0.0;
  for(int l=0; l<sz-1; l++) {
    h = x[l+1] - x[l];
    p = p + 0.5*(y[l] + y[l+1])*h;
  }

  return p;

}


void Trapsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y, Eigen::VectorXd& z)
{
  int sz = x.size();

  z[0] = 0.0;
  for(int l=1; l<sz; l++) {
    double h = x[l] - x[l-1];
    z[l] = z[l-1] + 0.5*(y[l-1] + y[l])*h;
  }
}

