// This may look like C code, but it is really -*- C++ -*-

// Interpolation routines operating on vectors

#ifndef _interp_h
#define _interp_h

#include <Vector.h>
// #include <Matrix.h>
#include <vector>

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

double odd2(double x, const Vector &xtab, const Vector &ftab, int even=0);

double odd2(double x, const vector<double> &xtab, const vector<double> &ftab, int even=0);

double drv2(double x, const Vector &xtab, const Vector &ftab, int even=0);

double drv2(double x, const vector<double> &xtab, const vector<double> &ftab, int even=0);


int Vlocate(double x, const Vector &xtab);

int Vlocate(double x, const vector<double> &xtab);


class Interp1d
{
public:

  virtual ~Interp1d();
  virtual double eval(const double& x);

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
  Spline1d &operator=(const Spline1d &);
  ~Spline1d();

  double eval(const double& x);
  void eval(const double& x, double& val, double& deriv);
  double deriv(const double& x);
};

class Linear1d : public Interp1d
{
private:

  Vector x;
  Vector y;

public:

  Linear1d();
  Linear1d(const Vector &x, const Vector &y);
  Linear1d &operator=(const Linear1d &);
  ~Linear1d();

  double eval(const double& x);
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
