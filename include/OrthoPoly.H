// Class for generating orhogonal polynomials

#ifndef _OrthoPoly_H
#define _OrthoPoly_H

#include <cfloat>
#include <cmath>
#include <Eigen/Eigen>

class OrthoPoly
{

 private:

  double coefs[4];

 protected:

  double a;
  double b;

  double *param;


 public:

  virtual double w(const double x) = 0;
  virtual double h(const int n) = 0;

  virtual double f0(const double x) = 0;
  virtual double f1(const double x) = 0;

  virtual double coef1(const int n) = 0;
  virtual double coef2(const int n) = 0;
  virtual double coef3(const int n) = 0;
  virtual double coef4(const int n) = 0;

  double f(const double x, const int n);
  Eigen::VectorXd fv(const double x, const int n);
  double get_a(void) { return a; }
  double get_b(void) { return b; }

};


class GenLagu : public OrthoPoly
{
  double alpha;

 public:
  
  GenLagu(double A) {
    a = 0.0;
    b = DBL_MAX;
    alpha = A; 
    param = &alpha;
  }

  double w(const double x) { return exp(-x)*pow(x, alpha); }
  double h(const int n) {return 1.0;}
  
  double coef1(const int n) { return 1.0+n; }
  double coef2(const int n) { return 1.0+alpha+2.0*n; }
  double coef3(const int n) { return -1.0; }
  double coef4(const int n) { return alpha + n; }

  inline double f0(const double x) { return 1.0; }
  inline double f1(const double x) { return -x + alpha + 1.0; }

};


class Legendre : public OrthoPoly
{

 public:
  
  Legendre(void) {
    a = -1.0;
    b =  1.0;
    param = 0;
  }

  double w(const double x) {return 1.0;}
  double h(const int n) { return 2.0/(2.0*n+1.0); }
  
  double coef1(const int n) { return 1.0+n; }
  double coef2(const int n) { return 0.0; }
  double coef3(const int n) { return 1.0+2.0*n; }
  double coef4(const int n) { return (double)n; }

  inline double f0(const double x) { return 1.0; }
  inline double f1(const double x) { return x; }

};


class Ultra : public OrthoPoly
{
  double alpha;

 public:
  
  Ultra(double A) {
    a = -1.0;
    b = 1.0;
    alpha = A;
    param = &alpha;
  }

  double w(const double x) { return pow(fabs(1.0 - x*x), alpha-0.5); }
  double h(const int n) 
    {if (alpha==0.0)
       return 2.0*M_PI/(n*n);
    else
      return M_PI*exp(log(2.0)*(1.0-2.0*alpha) + lgamma(2.0*alpha+n)
		      - lgamma(1.0+n) - 2.0*lgamma(alpha))/(alpha+n);
   }
  
  double coef1(const int n) { return 1.0+n; }
  double coef2(const int n) { return 0.0; }
  double coef3(const int n) { return 2.0*(alpha+n); }
  double coef4(const int n) { return 2.0*alpha + n - 1.0; }

  inline double f0(const double x) { return 1.0; }
  inline double f1(const double x) {
    if (alpha==0.0)
      return 2.0*x;
    else
      return 2.0*x*alpha;
  }

};


class Cheb1 : public OrthoPoly
{
 public:
  
  Cheb1(void) {
    a = -1.0;
    b = 1.0;
    param = 0;
  }
  
  double w(const double x) { return 1.0/sqrt(fabs(1.0 - x*x)); }
  double h(const int n) {if (n==0) return M_PI; else return 0.5*M_PI; }
    
  double coef1(const int n) { return 1.0; }
  double coef2(const int n) { return 0.0; }
  double coef3(const int n) { return 2.0; }
  double coef4(const int n) { return 1.0; }

  inline double f0(const double x) { return 1.0; }
  inline double f1(const double x) { return x; }

};

class Cheb2 : public OrthoPoly
{
 public:
  
  Cheb2(void) {
    a = -1.0;
    b = 1.0;
    param = 0;
  }
  
  double w(const double x) {return sqrt(fabs(1.0 - x*x)); }
  double h(const int n) {return 0.5*M_PI; }
    
  double coef1(const int n) { return 1.0; }
  double coef2(const int n) { return 0.0; }
  double coef3(const int n) { return 2.0; }
  double coef4(const int n) { return 1.0; }

  inline double f0(const double x) { return 1.0; }
  inline double f1(const double x) { return 2.0*x; }

};

class Hermite : public OrthoPoly
{
 public:
  
  Hermite(void) {
    a = -DBL_MAX;
    b =  DBL_MAX;
    param = 0;
  }
  
  double w(const double x) {return exp(-x*x); }
  double h(const int n) {return sqrt(M_PI)*exp(log(2.0)*n + lgamma(1.0+n)); }
    
  double coef1(const int n) { return 1.0; }
  double coef2(const int n) { return 0.0; }
  double coef3(const int n) { return 2.0; }
  double coef4(const int n) { return 2.0*n; }

  inline double f0(const double x) { return 1.0; }
  inline double f1(const double x) { return 2.0*x; }

};

#endif

