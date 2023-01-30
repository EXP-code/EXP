// Biorthnormal function class defintion for 1-d models
//

#ifndef _biorth1d_H
#define _biorth1d_H

#include <cmath>
#include <biorth.H>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

//! One-dimensional biorthgonal basis
class OneDBiorth : public Biorth
{
private:
  int dof;

public:
  //! Vertical wave number
  double kx;

  //! Constructor
  OneDBiorth(void) {dof = 1;}

  //! Get potential basis value
  virtual double potl(const int, const int, const double) = 0;

  //! Get density basis value
  virtual double dens(const int, const int, const double) = 0;

				// Fast for expansion! 
				// (using recursion relation)
  //! Write potential basis into a vetor
  virtual void potl(const int, const int, const double, Eigen::VectorXd&) = 0;
  //! Write density basis into a vetor
  virtual void dens(const int, const int, const double, Eigen::VectorXd&) = 0;

  //! Non-dimensional to dimensional coordinate map
  virtual double rb_to_r(double const) = 0;
  //! Dimensional to non-dimensional coordinate map
  virtual double r_to_rb(double const) = 0;
  //! Jacobian of mapping
  virtual double d_r_to_rb(double const) = 0;

  //! Minimum value in non-dimensional mapping
  virtual double rb_min(void) = 0;
  //! Maximum value in non-dimensional mapping
  virtual double rb_max(void) = 0;

  //! Number of degrees of freedom
  int get_dof(void) { return dof; }

  //! Return value of poential for given coefficient vector
  virtual double get_potl(double r, int l, Eigen::VectorXd& coef) = 0;

  //! Return value of density for given coefficient vector
  virtual double get_dens(double r, int l, Eigen::VectorXd& coef) = 0;

  //! Reset vertical wave number
  void reset_kx(double KX) { kx = KX; }
};

//! Basis based on trigonometric functions in vertical direction
class OneDTrig : public OneDBiorth
{
private:
  double zmax;
  int nrmax;
  Eigen::VectorXd kstar;
  Eigen::VectorXd kbstar;
  Eigen::VectorXd cnorm;
  void compute_kstar(int n);
  void compute_norm(void);

public:
  //! Enpoint offset for root finding (small)
  static double KSTOL;

  //! Root finding tolerance (small)
  static double KSZTOL;
  
  //! Constructor
  OneDTrig(void);

  //! Constructor
  OneDTrig(double kx);

  //! Constructor
  OneDTrig(double kx, double ZMAX);

  //! Reset grid with new parameters
  void reset(double KX, double ZMAX);

  //! Get potential value for basis
  double potl(const int n, const int tmp, const double z);

  //! Get density value for basis
  double dens(const int n, const int tmp, const double z);

  //! Get force value for basis
  double force(const int n, const int tmp, const double z);

				// Fast for expansion! 
				// (using recursion relation)
  //! Return potential basis in a vector
  void potl(const int n, const int tmp, const double z, Eigen::VectorXd& vec);

  //! Return density basis in a vector
  void dens(const int n, const int tmp, const double z, Eigen::VectorXd& vec);

  //! Return force basis in a vector
  void force(const int n, const int tmp, const double z, Eigen::VectorXd& vec);

  //! Non-dimensional to dimensional mapping
  double rb_to_r(double const x)    { return x; }

  //! Dimensional to non-dimensional mapping
  double r_to_rb(double const x)    { return x; }

  //! Jacobian for non-dimensional mapping
  double d_r_to_rb(double const x)  { return 1.0; }

  //! Minimum value in non-dimensional mapping
  double rb_min(void)               { return -zmax; }

  //! Maximum value in non-dimensional mapping
  double rb_max(void)               { return zmax; }

  //! Kernel value
  double krnl(int n, int k=0)         { return  1.0; }

  //! Norm value
  double norm(int n, int k=0)         { return  1.0; }

  //! Return potential field for given coefficient vector
  double get_potl(double r, int l, Eigen::VectorXd& coef);

  //! Return density field for given coefficient vector
  double get_dens(double r, int l, Eigen::VectorXd& coef);

  //! Return force field for given coefficient vector
  double get_force(double r, int l, Eigen::VectorXd& coef);

};

#endif


