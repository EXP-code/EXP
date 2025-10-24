/*
 *	Computes integral using spline fit
 *	with a supplied integrand
 *
 *	Synopsis:
 *		ans = splsum(x,y,n);
 *		double	x[],y[];
 *              int     n;
 *
 * 	Use:
 * 		x       domain
 * 		y       range
 * 		n       number of values y[1]-->y[n]
 * 
 *      Notes:
 *              If (n<4), trapezoidal rule is used
 *
 */

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "interp.H"

double	Splsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y)
{
  static Eigen::VectorXd y2;

  int n = x.size();
  
  if (n < 2) {
    std::cerr << "Splsum error, can't do intgral with one grid point!"
	      << std::endl << std::endl;
    exit (-1);
  }
  else if (n < 3)
    return 0.5*( (x[1]-x[0])*(y[0]+y[1]) );
  else if (n < 4)
    return 0.5*( (x[1]-x[0])*(y[0]+y[1]) + 
		 (x[2]-x[1])*(y[1]+y[2]) );

  y2.resize(n);
  Spline(x, y, -1.0e30, -1.0e30, y2);

  double h, p=0.0;
  for(int l=0; l<n-1; l++) {
    h = x[l+1] - x[l];
    p = p + 0.5*(y[l] + y[l+1])*h - (y2[l] + y2[l+1])*h*h*h/24.0;
  }

  return p;
}


void Splsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y, Eigen::VectorXd& z)
{
  static Eigen::VectorXd y2;

  int n = x.size();

  if (n < 2) {
    std::cerr << "Splsum error, can't do intgral with one grid point!"
	      << std::endl << std::endl;
    z = Eigen::VectorXd::Zero(n);
    return;
  }
  else if (n < 4) {
    z[0] = 0.0;
    for (int l=1; l<n; l++)
      z[l] = z[l-1] + 0.5*(y[l-1]+y[l])*(x[l]-x[l-1]);
    return;
  }

  y2.resize(n);
  Spline(x, y, -1.0e30, -1.0e30, y2);

  double h;

  z[0] = 0.0;
  for(int l=1; l<n; l++) {
    h = x[l] - x[l-1];
    z[l] = z[l-1] + 0.5*(y[l-1] + y[l])*h - (y2[l-1] + y2[l])*h*h*h/24.0;
  }
}

