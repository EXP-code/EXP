// This may look like C code, but it is really -*- C++ -*-

// Biorthnormal function class defintion


#ifndef _biorth2d_H
#define _biorth2d_H

#include <string>
#include "biorth.H"

class AxiSymBiorth;

class CBDisk : public AxiSymBiorth
{
private:
  int numz;
  double phif(const int n, const int m, const double r);

public:

  CBDisk(void);

				// Required functions

  double potl(const int n, const int m, const double r);
  double dens(const int n, const int m, const double r);
  inline double krnl(const int n, const int m) {return 1.0;}
  double norm(const int n, const int m);
  void potl(const int n, const int m, const double r, Eigen::VectorXd& a);
  void dens(const int n, const int m, const double r, Eigen::VectorXd& a);


				// Supplemental functions

  inline double potlR(const int n, const int m, const double r) {
    return potl(n, m, r); }
  inline double densR(const int n, const int m, const double r) {
    return dens(n, m, r); }

  double potlRZ(const int n, const int m, const double r, const double z);
  void set_numz(int n) { numz = n;}

  inline double rb_to_r(double const r) {return r;}
  inline double r_to_rb(double const r) {return r;}
  double d_r_to_rb(double const r) {return 1.0;}

  double rb_min(void) { return 0.0; }
  double rb_max(void) { return 1.0e30; }

};




class BSDisk : public AxiSymBiorth
{
private:
  int numz;
  int nmax, mmax;
  double rmax;

  std::vector<Eigen::VectorXd> a;

  Eigen::MatrixXd t_f, t_y;
  Eigen::VectorXd t_dr;
  int t_n;

  void setup_potl_table(void);

public:

  BSDisk(double RMAX, int NMAX, int MMAX);

  ~BSDisk();
				// Required functions

  double potl(const int n, const int m, const double r);
  double dens(const int n, const int m, const double r);
  inline double krnl(const int n, const int m) {return 1.0;}
  double norm(const int n, const int m);
  void potl(const int n, const int m, const double r, Eigen::VectorXd& a);
  void dens(const int n, const int m, const double r, Eigen::VectorXd& a);


				// Supplemental functions

  inline double potlR(const int n, const int m, const double r) {
    return potl(n, m, r); }
  inline double densR(const int n, const int m, const double r) {
    return dens(n, m, r); }

  double potlRZ(const int n, const int m, const double r, const double z);

  inline double rb_to_r(double const r) {return r;}
  inline double r_to_rb(double const r) {return r;}
  double d_r_to_rb(double const r) {return 1.0;}

  double rb_min(void) { return 0.0; }
  double rb_max(void) { return rmax; }

};

enum BiorthFcts2d {clutton_brock_2d, bessel_2d};
// static const char* BiorthFcts2dName[] = {"CBDisk", "BSDisk"};

#endif
