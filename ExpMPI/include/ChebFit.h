// -*- C++ -*-

//
// Chebyshev fitting and smoothing class
//


#ifndef _ChebFit_h

#define _ChebFit_h 1

class Vector;

class ChebFit
{
private:
  int n;
  double a, b;
  Vector c;
  Vector c1;
  Vector c2;
  
  void chder(Vector& cin, Vector& cder);
  double chebev(double x, Vector& cin);

  bool rcheck;
  bool defined;

public:
				// Constructors
  
  ChebFit(void) : rcheck(true), defined(false) {};
  ChebFit(double A, double B, Vector& X, Vector &Y, int N);
  ChebFit(Vector &X, Vector &Y, int N);

				// Member functions

  void range_check_on() { rcheck = true; }
  void range_check_off() { rcheck = false; }
  void new_limits(double A, double B);
  void new_func(double (*func)(double), double A, double B, int N);
  void new_data(Vector &X, Vector &Y, int N);

  double eval(double x) 
  {if (defined) return chebev(x, c ); else bomb("no data!"); return 0.0;}
  double der1(double x) 
  {if (defined) return chebev(x, c1); else bomb("no data!"); return 0.0;}
  double der2(double x) 
  {if (defined) return chebev(x, c2); else bomb("no data!"); return 0.0;}
  
  void bomb(char *, ...);
  
};

#endif
