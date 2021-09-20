//
// Chebyshev fitting and smoothing class
//

#ifndef _ChebFit_H
#define _ChebFit_H

#include <Eigen/Eigen>

class ChebFit
{
private:
  int n;
  double a, b;
  Eigen::VectorXd c;
  Eigen::VectorXd c1;
  Eigen::VectorXd c2;
  
  void chder(Eigen::VectorXd& cin, Eigen::VectorXd& cder);
  double chebev(double x, Eigen::VectorXd& cin);

  bool rcheck;
  bool defined;

public:
				// Constructors
  
  ChebFit(void) : rcheck(true), defined(false) {};
  ChebFit(double A, double B, Eigen::VectorXd& X, Eigen::VectorXd &Y, int N);
  ChebFit(Eigen::VectorXd &X, Eigen::VectorXd &Y, int N);

				// Member functions

  void range_check_on() { rcheck = true; }
  void range_check_off() { rcheck = false; }
  void new_limits(double A, double B);
  void new_func(double (*func)(double), double A, double B, int N);
  void new_data(Eigen::VectorXd &X, Eigen::VectorXd &Y, int N);

  double eval(double x) 
  {if (defined) return chebev(x, c ); else bomb("no data!"); return 0.0;}
  double der1(double x) 
  {if (defined) return chebev(x, c1); else bomb("no data!"); return 0.0;}
  double der2(double x) 
  {if (defined) return chebev(x, c2); else bomb("no data!"); return 0.0;}
  
  void bomb(const char *, ...);
  
};

#endif
