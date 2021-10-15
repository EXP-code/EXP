#ifndef _interp_H
#define _interp_H

#include <Eigen/Eigen>
#include <vector>
#include <deque>

void Spline(const Eigen::VectorXd &x, const Eigen::VectorXd &y, double yp1, double ypn, Eigen::VectorXd &y2);

void Splint1(const Eigen::VectorXd &xa, const Eigen::VectorXd &ya, const Eigen::VectorXd &y2a, double x, 
	     double &y, int even=0);

void Splint2(const Eigen::VectorXd &xa, const Eigen::VectorXd &ya, const Eigen::VectorXd &y2a, double x, 
	     double &y, double &dy, int even=0);

void Splint3(const Eigen::VectorXd &xa, const Eigen::VectorXd &ya, const Eigen::VectorXd &y2a, double x, 
	     double &y, double &dy, double &dyy, int even=0);

double Splsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y);
void Splsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y, Eigen::VectorXd& z);

double Trapsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y);
void Trapsum(const Eigen::VectorXd& x, const Eigen::VectorXd& y, Eigen::VectorXd& z);

// This will work for both STL deque and vector
//
template <class V>
double odd2(double x, const V& xtab, const V& ftab, int even=0);

template <class V>
double drv2(double x, const V& xtab, const V& ftab, int even=0);

template <class V>
int Vlocate(double x, const V& xtab);

template <class V>
int Vlocate_with_guard(double value, const V& vec);

class Interp1d
{
public:

  virtual ~Interp1d() = 0;
  virtual double eval(const double& x) = 0;
  virtual double deriv(const double& x) = 0;

};

class Spline1d : public Interp1d
{
private:

  Eigen::VectorXd x;
  Eigen::VectorXd y, y2;

public:

  Spline1d();
  Spline1d(const Eigen::VectorXd &x, const Eigen::VectorXd &y, 
	   double d1=-1.0e30, double d2=-1.0e30);
  Spline1d(const std::vector<double> &x, const std::vector<double> &y, 
	   double d1=-1.0e30, double d2=-1.0e30);
  Spline1d &operator=(const Spline1d &);
  ~Spline1d();

  double eval(const double& x);
  void eval(const double& x, double& val, double& deriv);
  double deriv(const double& x);
};


class Cheby1d : public Interp1d
{
private:

  int n;
  double a, b;
  std::vector<double> c;
  std::vector<double> c1;
  std::vector<double> c2;
  
  void   chder(std::vector<double>& cin, std::vector<double>& cder);
  double chebev(double x, std::vector<double>& cin);
  double chebint(double x);

  bool defined;


public:

  Cheby1d();
  Cheby1d(const double A, const double B, 
	  const Eigen::VectorXd &x, const Eigen::VectorXd &y, int N);
  Cheby1d(const double A, const double B, 
	  std::vector<double> &x, std::vector<double> &y, int N);
  Cheby1d(const Eigen::VectorXd &x, const Eigen::VectorXd &y, int N);
  Cheby1d(std::vector<double> &x, std::vector<double> &y, int N);
  Cheby1d &operator=(const Cheby1d &);
  ~Cheby1d();

  double eval(const double& x)
  {if (defined) return chebev(x, c);  else bomb("no data!"); return 0.0;}
  double deriv(const double& x)
  {if (defined) return chebev(x, c1); else bomb("no data!"); return 0.0;}
  double integral(const double& x)
  {if (defined) return chebint(x);    else bomb("no data!"); return 0.0;}

				// New member functions

  void new_limits(double A, double B);
  void new_func(double (*func)(double), double A, double B, int N);
  void new_data(std::vector<double> &X, std::vector<double> &Y, int N);

  void eval(const double& x, double& val, double& deriv) {
    if (defined) {
      val = chebev(x, c); 
      deriv = chebev(x, c1);
    } else {
      val = deriv = 0.0;
      bomb("no data!");
    }
  }

  double deriv2(const double& x)
  {if (defined) return chebev(x, c1); else bomb("no data!"); return 0.0;}
  void bomb(const char *, ...);
};


class Linear1d : public Interp1d
{
private:

  Eigen::VectorXd x;
  Eigen::VectorXd y;

public:

  Linear1d();
  Linear1d(const Eigen::VectorXd &x, const Eigen::VectorXd &y);
  Linear1d(const std::vector<double> &x, const std::vector<double> &y);
  Linear1d &operator=(const Linear1d &);
  ~Linear1d();

  double eval(const double& x);
  double deriv(const double& x);
};


class Spline2d
{
private:

  Eigen::VectorXd x, y;
  Eigen::VectorXd xval, xval2;
  Eigen::MatrixXd mat, mat2d;

public:
  static double DERIV;

  Spline2d(void);
  Spline2d(const Eigen::VectorXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &mat);
  Spline2d &operator=(const Spline2d &);

  double eval(const double& x, const double& y);
};

#include <MonotCubicInterpolator.H>

#endif
