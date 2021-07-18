// This may look like C code, but it is really -*- C++ -*-

// Interpolation routines operating on vectors

#ifndef _interp_h
#define _interp_h

#include <Vector.h>
#include <vector>
#include <deque>

void Spline(const Vector &x, const Vector &y, double yp1, double ypn, Vector &y2);

void Splint1(const Vector &xa, const Vector &ya, const Vector &y2a, double x, 
	     double &y, int even=0);

void Splint2(const Vector &xa, const Vector &ya, const Vector &y2a, double x, 
	     double &y, double &dy, int even=0);

void Splint3(const Vector &xa, const Vector &ya, const Vector &y2a, double x, 
	     double &y, double &dy, double &dyy, int even=0);

double Splsum(const Vector& x, const Vector& y);
void Splsum(const Vector& x, const Vector& y, Vector& z);

double Trapsum(const Vector& x, const Vector& y);
void Trapsum(const Vector& x, const Vector& y, Vector& z);

// For old Vector class
//
double odd2(double x, const Vector & xtab, const Vector &ftab, int even=0);

// This will work for both STL deque and vector
//
template <class V>
double odd2(double x, const V& xtab, const V& ftab, int even=0);

double drv2(double x, const Vector& xtab, const Vector &ftab, int even=0);
template <class V>
double drv2(double x, const V& xtab, const V& ftab, int even=0);

int Vlocate(double x, const Vector& xtab);

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

  Vector x;
  Vector y, y2;

public:

  Spline1d();
  Spline1d(const Vector &x, const Vector &y, 
	   double d1=-1.0e30, double d2=-1.0e30);
  Spline1d(const vector<double> &x, const vector<double> &y, 
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
  vector<double> c;
  vector<double> c1;
  vector<double> c2;
  
  void chder(vector<double>& cin, vector<double>& cder);
  double chebev(double x, vector<double>& cin);

  bool defined;


public:

  Cheby1d();
  Cheby1d(const double A, const double B, 
	  const Vector &x, const Vector &y, int N);
  Cheby1d(const double A, const double B, 
	  vector<double> &x, vector<double> &y, int N);
  Cheby1d(const Vector &x, const Vector &y, int N);
  Cheby1d(vector<double> &x, vector<double> &y, int N);
  Cheby1d &operator=(const Cheby1d &);
  ~Cheby1d();

  double eval(const double& x)
  {if (defined) return chebev(x, c); else bomb("no data!"); return 0.0;}
  double deriv(const double& x)
  {if (defined) return chebev(x, c1); else bomb("no data!"); return 0.0;}

				// New member functions

  void new_limits(double A, double B);
  void new_func(double (*func)(double), double A, double B, int N);
  void new_data(vector<double> &X, vector<double> &Y, int N);

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

  Vector x;
  Vector y;

public:

  Linear1d();
  Linear1d(const Vector &x, const Vector &y);
  Linear1d(const vector<double> &x, const vector<double> &y);
  Linear1d &operator=(const Linear1d &);
  ~Linear1d();

  double eval(const double& x);
  double deriv(const double& x);
};


class Spline2d
{
private:

  Vector x, y;
  Vector xval, xval2;
  Matrix mat, mat2d;

public:
  static double DERIV;

  Spline2d(void);
  Spline2d(const Vector &x, const Vector &y, const Matrix &mat);
  Spline2d &operator=(const Spline2d &);

  double eval(const double& x, const double& y);
};

#endif
