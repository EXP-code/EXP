#include <cmath>
#include <vector>

#include "interp.H"

Interp1d::~Interp1d(void) 
{
				// Nothing
}

double Interp1d::eval(const double& x) 
{
				// Nothing
  return 0.0;
}

Linear1d::Linear1d()
{
				// Nothing
}

Linear1d::~Linear1d()
{
				// Nothing
}

Linear1d &Linear1d::operator=(const Linear1d &p)
{
  x = p.x;
  y = p.y;

  return *this;
}

Linear1d::Linear1d(const Eigen::VectorXd &X, const Eigen::VectorXd &Y)
{
  x = X;
  y = Y;
}


Linear1d::Linear1d(const std::vector<double> &X, const std::vector<double> &Y)
{
  int sz = X.size();
  x.resize(sz);
  y.resize(sz);

  for (int i=0; i<sz; i++) {
    x[i] = X[i];
    y[i] = Y[i];
  }
}


double Linear1d::eval(const double &x1)
{
  return odd2(x1, x, y);
}

double Linear1d::deriv(const double& x1)
{
  return drv2(x1, x, y);
}


Spline1d::Spline1d()
{
				// Nothing
}

Spline1d::~Spline1d() 
{
				// Nothing
}


Spline1d &Spline1d::operator=(const Spline1d &p)
{
  x = p.x;
  y = p.y;
  y2 = p.y2;

  return *this;
}

Spline1d::Spline1d(const Eigen::VectorXd &X, const Eigen::VectorXd &Y, double d1, double d2)
{
  x = X;
  y = Y;
  y2.resize(Y.size());

  Spline(x, y, d1, d2, y2);
}


Spline1d::Spline1d(const std::vector<double> &X, const std::vector<double> &Y, 
		   double d1, double d2)
{
  int sz = X.size();
  x .resize(sz);
  y .resize(sz);
  y2.resize(sz);

  for (int i=0; i<sz; i++) {
    x[i] = X[i];
    y[i] = Y[i];
  }

  Spline(x, y, d1, d2, y2);
}


double Spline1d::eval(const double &x1)
{
  double ans;
  Splint1(x, y, y2, x1, ans);
  return ans;
}

void Spline1d::eval(const double &x1, double& val, double &deriv)
{
  Splint2(x, y, y2, x1, val, deriv);
}

double Spline1d::deriv(const double &x1)
{
  double ans, dum;
  Splint2(x, y, y2, x1, dum, ans);
  return ans;
}
