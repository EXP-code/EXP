// This may look like C code, but it is really -*- C++ -*-

// Biorthnormal function class defintion


#ifndef _biorth1d_h
#ifdef __GNUG__
#pragma interface
#endif

#define _biorth1d_h 1

static const char rcsid_biorth1d[] = "$Id$";

#include <math.h>
#include <biorth.h>

class OneDBiorth : public Biorth
{
private:
  int dof;

public:
  double kx;

  OneDBiorth(void) {dof = 1;}

  virtual double potl(const int, const int, const double) = 0;
  virtual double dens(const int, const int, const double) = 0;

				// Fast for expansion! 
				// (using recursion relation)
  virtual void potl(const int, const int, const double, Vector&) = 0;
  virtual void dens(const int, const int, const double, Vector&) = 0;

  virtual double rb_to_r(double const) = 0;
  virtual double r_to_rb(double const) = 0;
  virtual double d_r_to_rb(double const) = 0;

  virtual double rb_min(void) = 0;
  virtual double rb_max(void) = 0;

  int get_dof(void) { return dof; }

  virtual double get_potl(double r, int l, Vector& coef) = 0;
  virtual double get_dens(double r, int l, Vector& coef) = 0;

  void reset_kx(double KX) { kx = KX; }
};

class OneDTrig : public OneDBiorth
{
private:
  double zmax;
  int nrmax;
  Vector kstar;
  Vector kbstar;
  Vector cnorm;
  void compute_kstar(int n);
  void compute_norm(void);

public:
  static double KSTOL;
  static double KSZTOL;
  
				// Constructors
  OneDTrig(void);
  OneDTrig(double kx);
  OneDTrig(double kx, double ZMAX);

  void reset(double KX, double ZMAX);

  double potl(const int n, const int tmp, const double z);
  double dens(const int n, const int tmp, const double z);
  double force(const int n, const int tmp, const double z);

				// Fast for expansion! 
				// (using recursion relation)
  void potl(const int n, const int tmp, const double z, Vector& vec);
  void dens(const int n, const int tmp, const double z, Vector& vec);
  void force(const int n, const int tmp, const double z, Vector& vec);

  double rb_to_r(double const x)    { return x; }
  double r_to_rb(double const x)    { return x; }
  double d_r_to_rb(double const x)  { return 0.0; }

  double rb_min(void)               { return 0.0; }
  double rb_max(void)               { return zmax; }

  double krnl(int n, int k=0)         { return  1.0; }
  double norm(int n, int k=0)         { return  1.0; }

  double get_potl(double r, int l, Vector& coef);
  double get_dens(double r, int l, Vector& coef);
  double get_force(double r, int l, Vector& coef);

};

#endif


