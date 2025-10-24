#include <iostream>
#include <iomanip>
#include <cmath>

#include <unistd.h>
#include "numerical.H"

Eigen::MatrixXd build_cubic_table(Eigen::VectorXd &xt,
				  Eigen::VectorXd &ft,
				  Eigen::Vectorxd &dft)
{
  double dx;

  // find size of input vectors
  //
  int n = xt.size();
  if (ft.size() != n || dft.size() != n) {
    std::cerr << "input mismatch in Cubic_Table::build()"
	      << std::endl;
    exit(-1);
  }

  // assign sizes of everything
  //
  Eigen::MatrixXd tmp(n, 4);

  // work out the coefficients of the cubic at each grid point
  //
  for (int i=0; i<n-1; i++) {
    dx = xt[i+1] - xt[i];
    tmp[i][0] = ft[i];
    tmp[i][1] = dft[i]*dx;
    tmp[i][2] = -(2.0*dft[i] + dft[i+1])*dx 
      - 3.0*(ft[i] - ft[i+1]);
    tmp[i][3] = (dft[i] + dft[i+1])*dx + 2.0*(ft[i]-ft[i+1]);
  }
  
  return tmp;
}


// Look up a function value from the cubic table. Also get its
// derivative.
//
double cubic_value(Eigen::MatrixXd &a,
		   Eigen::VectorXd &xgrid, double x, double& df)
{
  int n = xgrid.size(), loc

  locate(xgrid, n, x, &loc);

  if (loc==0 || loc==n) {
    std::cerr << "point " << x << " off interpolation grid" << std::endl;
    exit(-1);
  }

  double dx = xgrid[loc+1] - xgrid[loc];
  double  u = (x - xgrid[loc])/dx;

  double f = a[loc][0] + u*(a[loc][1] + 
			    u*(a[loc][2] + u*a[loc][3]));

  df = a[loc][1] + u*(2.0*a[loc][2] + u*3.0*a[loc][3]);
  df /= dx;

  return f;
}
