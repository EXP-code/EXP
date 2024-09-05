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

  //! Destructor
  virtual ~Interp1d() = 0;

  //! Evaluate the interpolant
  virtual double eval(const double& x) = 0;

  //! Evaluate first derivative
  virtual double deriv(const double& x) = 0;

  //@{
  //! Data limits
  virtual double xlo() = 0;
  virtual double xhi() = 0;
  //@}
};

class Spline1d : public Interp1d
{
private:

  Eigen::VectorXd x;
  Eigen::VectorXd y, y2;

public:

  //! Null constructor
  Spline1d();

  //! Consruct with Eigen input. d1 and d2 are 2nd derivative BCs
  Spline1d(const Eigen::VectorXd &x, const Eigen::VectorXd &y, 
	   double d1=-1.0e30, double d2=-1.0e30);

  //! Consruct with std::vector input. d1 and d2 are 2nd derivative BCs
  Spline1d(const std::vector<double> &x, const std::vector<double> &y, 
	   double d1=-1.0e30, double d2=-1.0e30);

  //! Copy constructor
  Spline1d &operator=(const Spline1d &);

  //! Destructor
  ~Spline1d();

  //! Evaluate the interpolant
  double eval(const double& x);

  //! Evaluate the interpolant and its first derivative
  void eval(const double& x, double& val, double& deriv);

  //! Evaluate first derivative
  double deriv(const double& x);

  //@{
  //! Data limits
  double xlo() { return x(0); }
  double xhi() { return x(x.size()-1); }
  //@}
};


class Cheby1d : public Interp1d
{
private:

  int n;
  double a, b, xmin=0, xmax=0;
  std::vector<double> c;
  std::vector<double> c1;
  std::vector<double> c2;
  
  void   chder(std::vector<double>& cin, std::vector<double>& cder);
  double chebev(double x, std::vector<double>& cin);
  double chebint(double x);

  bool defined;


public:

  //! Null constructor
  Cheby1d();

  //! Construct with Eigen input (N is Chebyshev order)
  Cheby1d(const double A, const double B, 
	  const Eigen::VectorXd &x, const Eigen::VectorXd &y, int N);

  //! Construct with std::vector input (N is Chebyshev order)
  Cheby1d(const double A, const double B, 
	  std::vector<double> &x, std::vector<double> &y, int N);

  //! Construct with Eigen input (N is Chebyshev order) using data limits
  Cheby1d(const Eigen::VectorXd &x, const Eigen::VectorXd &y, int N);

  //! Construct with std::vector  input (N is Chebyshev order) using data limits
  Cheby1d(std::vector<double> &x, std::vector<double> &y, int N);

  //! Copy constructor
  Cheby1d &operator=(const Cheby1d &);

  //! Destructor
  ~Cheby1d();

  //! Evaluate the interpolant
  double eval(const double& x)
  {if (defined) return chebev(x, c);  else bomb("no data!"); return 0.0;}

  //! Evaluate first derivative
  double deriv(const double& x)
  {if (defined) return chebev(x, c1); else bomb("no data!"); return 0.0;}

  //! Evaluate integral of interpolant
  double integral(const double& x)
  {if (defined) return chebint(x);    else bomb("no data!"); return 0.0;}

				// New member functions

  //! Update the data limits
  void new_limits(double A, double B);

  //! Generate interpolant from a function
  void new_func(std::function<double(double)> func, double A, double B, int N);

  //! Update the data
  void new_data(std::vector<double> &X, std::vector<double> &Y, int N);

  //! Evaluate the interpolant and its first derivative
  void eval(const double& x, double& val, double& deriv) {
    if (defined) {
      val = chebev(x, c); 
      deriv = chebev(x, c1);
    } else {
      val = deriv = 0.0;
      bomb("no data!");
    }
  }

  //! Evaluate the second derivative
  double deriv2(const double& x)
  {if (defined) return chebev(x, c1); else bomb("no data!"); return 0.0;}

  //! Handle error
  void bomb(const char *, ...);

  //@{
  //! Data limits
  double xlo() { return xmin; }
  double xhi() { return xmax; }
  //@}
};


class Linear1d : public Interp1d
{
private:

  Eigen::VectorXd x;
  Eigen::VectorXd y;

public:

  //! Null constructor
  Linear1d();

  //! Construct from Eigen input
  Linear1d(const Eigen::VectorXd &x, const Eigen::VectorXd &y);

  //! Construct from std::vector input
  Linear1d(const std::vector<double> &x, const std::vector<double> &y);

  //! Copy constructor
  Linear1d &operator=(const Linear1d &);

  //! Destructor
  ~Linear1d();

  //! Evaluate the interpolant
  double eval(const double& x);

  //! Evaluate first derivative
  double deriv(const double& x);

  //@{
  //! Data limits
  double xlo() { return x(0); }
  double xhi() { return x(x.size()-1); }
  //@}
};


//! Two dimensional spline interpolator
class Spline2d
{
private:

  Eigen::VectorXd x, y;
  Eigen::VectorXd xval, xval2;
  Eigen::MatrixXd mat, mat2d;

public:
  static double DERIV;

  //! Null constructor
  Spline2d(void);

  //! Construct from Eigen input
  Spline2d(const Eigen::VectorXd &x, const Eigen::VectorXd &y, const Eigen::MatrixXd &mat);

  //! Copy constructor
  Spline2d &operator=(const Spline2d &);

  //! Evaluate the interpolant
  double eval(const double& x, const double& y);
};

#include <MonotCubicInterpolator.H>

#endif
