#ifndef _sphereSL_h
#define _sphereSL_h

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include <numerical.H>
#include <biorth.H>

#include <SLGridMP2.H>


class SphereSL : public AxiSymBiorth
{

private:
  
  int lmax, nmax, numr;
  double rmin, rmax, scale;

  SLGridSph *slgrid;
  

public:

				// Global MPI stuff
  static int mpi;		// default: 0=off

				// Check for cached table
  static bool cache;		// default: true

				// Constructors

  SphereSL(int lmax, int nmax, int numr, double rmin, double rmax,
	   double scale, std::shared_ptr<SphericalModelTable> mod);
  ~SphereSL(void);

				// Required members

  double potl(const int n, const int l, const double x);
  double dens(const int n, const int l, const double x);

  inline double potlR(const int n, const int l, const double r) {
    return potl(n, l, r_to_rb(r) ); }

  inline double potlRZ(const int n, const int l, const double r,
		       const double z) {
    return potl(n, l, r_to_rb(sqrt(r*r+z*z)) ); }


  void potl(const int n, const int l, const double x, Eigen::VectorXd& t);
  void dens(const int n, const int l, const double x, Eigen::VectorXd& t);

  

  double r_to_rb(double const r);
  double rb_to_r(double const x);
  double d_r_to_rb(double const x);
  double rb_min(void) { return slgrid->r_to_xi(rmin); } 
  double rb_max(void) { return slgrid->r_to_xi(rmax); }

  double get_potl(const double r, const int l, const Eigen::VectorXd& coef);
  double get_dens(const double r, const int l, const Eigen::VectorXd& coef);

  inline double krnl(const int n, const int l) { return 1.0; }
  inline double norm(const int n, const int l) { return 1.0; }

  double d_rb_to_r(const double x);

};

#endif // _BiorthSL_h

