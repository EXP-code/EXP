//
// Chebyshev fitting and smoothing class
//

#include <functional>
#include <cstdarg>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

#include <interp.H>


Cheby1d::Cheby1d()
{
  defined = false;
}


Cheby1d::~Cheby1d() 
{
				// Nothing
}

Cheby1d &Cheby1d::operator=(const Cheby1d &p)
{
  n = p.n;
  a = p.a;
  b = p.b;
  c = p.c;
  c1 = p.c1;
  c2 = p.c2;
  defined = p.defined;
  return *this;
}


Cheby1d::Cheby1d(std::vector<double>& X, std::vector<double> &Y, int N)
{

				// set private values
  n = N;
  a = X[0];
  b = X[X.size()-1];

  new_data(X, Y, N);

}

Cheby1d::Cheby1d(double A, double B, 
		 std::vector<double>& X, std::vector<double> &Y, int N)
{
  
				// set private values
  n = N;
  a = A;
  b = B;

  new_data(X, Y, N);

}


Cheby1d::Cheby1d(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, int N)
{

				// set private values
  n = N;
  a = X[0];
  b = X[X.size()-1];

  std::vector<double> XX(X.size()), YY(X.size());
  for (int n=0; n<X.size(); n++) {
    XX[n] = X[n];
    YY[n] = Y[n];
  }

  new_data(XX, YY, N);

}

Cheby1d::Cheby1d(double A, double B, 
		 const Eigen::VectorXd& X, const Eigen::VectorXd& Y, int N)
{
  
				// set private values
  n = N;
  a = A;
  b = B;

  std::vector<double> XX(X.size()), YY(X.size());
  for (int n=0; n<X.size(); n++) {
    XX[n] = X[n];
    YY[n] = Y[n];
  }

  new_data(XX, YY, N);
}


void Cheby1d::new_limits(double A, double B)
{
  a = A;
  b = B;
}
    

void Cheby1d::new_func(std::function<double(double)> func,
		       double A, double B, int N)
{
  double y, sum, fac, bpa, bma;

  a = A;
  b = B;
  n = N;

  std::vector<double> f(n);
  c  = std::vector<double>(n);
  c1 = std::vector<double>(n);
  c2 = std::vector<double>(n);

  bma = 0.5*(b - a);
  bpa = 0.5*(b + a);
  for (int k=0; k<n; k++) {
    y = bma*cos(M_PI*(k+0.5)/n) + bpa;
    f[k] = func(y);
  }

  fac=2.0/n;
  for (int j=0; j<n; j++) {
    sum=0.0;
    for (int k=0;k<n;k++)
      sum += f[k]*cos(M_PI*j*(k+0.5)/n);
    c[j]=fac*sum;
  }

  chder(c, c1);
  chder(c1, c2);

  defined = true;
}


void Cheby1d::new_data(std::vector<double>& X, std::vector<double>& Y, int N)
{
  double y, sum, fac, bpa, bma;

  std::vector<double> f(n);

  c  = std::vector<double>(n);
  c1 = std::vector<double>(n);
  c2 = std::vector<double>(n);

  bma = 0.5*(b - a);
  bpa = 0.5*(b + a);
  for (int k=0; k<n; k++) {
    y = bma*cos(M_PI*(k+0.5)/n) + bpa;
    f[k] = odd2(y, X, Y);
  }

  fac=2.0/n;
  for (int j=0; j<n; j++) {
    sum=0.0;
    for (int k=0; k<n; k++)
      sum += f[k]*cos(M_PI*j*(k+0.5)/n);
    c[j]=fac*sum;
  }

  chder(c,  c1);
  chder(c1, c2);

  defined = true;
}


void Cheby1d::chder(std::vector<double>& cin, std::vector<double>& cder)
{
  double con;

  cder[n-1] = 0.0;
  cder[n-2] = 2*(n-1)*cin[n-1];
  for (int j=n-3; j>=0; j--)
    cder[j] = cder[j+2]+ 2*(j+1)*cin[j+1];

  con = 2.0/(b-a);
  for (int j=0; j<n; j++)
    cder[j] *= con;
}


double Cheby1d::chebev(double x, std::vector<double>& cin)
{
  double d=0.0, dd=0.0, sv, y, y2;
  if ((x-a)*(x-b) > 0.0) bomb("Cheby1d::chebev", "x out of range", 0);

  y2 = 2.0*(y=(2.0*x-a-b)/(b-a));

  for (int j=n-1;j>=1;j--) {
    sv = d;
    d = y2*d-dd+cin[j];
    dd = sv;
  }

  return y*d-dd+0.5*cin[0];
}

double Cheby1d::chebint(double x)
{
  if ((x-a)*(x-b) > 0.0) bomb("Cheby1d::integral", "x out of range", 0);

  double y = (2.0*x-a-b)/(b-a);

  // Evaluate by recursion
  //
  std::vector<double> t(n+1);
  t[0] = 1.0;
  t[1] = y;
  for (int j=2; j<n+1; j++) t[j] = 2.0*t[j-1] - t[j-2];

  // Integral term by term
  //
				// First two terms
  double ret = c[0]*(y + 1.0) + c[1]*0.5*(y*y - 1.0);
  double tm1 = -1.0;		// End-point eval
  for (int j=2; j<n; j++) {	// Remaining terms
    ret += (t[j+1] - tm1)/(j + 1) - (t[j-1] - tm1)/(j-1);
    tm1 *= -1.0;		// Update end-point eval
  }  

  return 0.5*(b - a)*ret;	// The answer . . . 
}

void Cheby1d::bomb(const char *a, ...)
{
  va_list ap;
  char *b, *c;

  va_start(ap, a);
  b = va_arg(ap, char *);
  c = va_arg(ap, char *);
  va_end(ap);

  if (b != (char *)0)
    std::cerr << a << ": " << b;
  if (c != (char *)0)
    std::cerr << ", " << c;
  std::cerr << std::endl;

  exit(-1);
}

