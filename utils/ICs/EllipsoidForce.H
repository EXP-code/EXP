#ifndef EllipForce_H
#define EllipForce_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>

#include <Eigen/Eigen>

#include <interp.H>
#include <gaussQ.H>

/** Defines potential and force for ellipsoid models
    
    @param a is the major axis
    @param b is the submajor (middle) axis
    @param c is the minor axis
    @param mass is the total mass for the ellipsoid
    @param num is the number of integration knots
    @param numr is the number of grid points in the radial table

    \note{Mathematics}

    Internal potential has the form (Chandra pg. 52, eqn. 93)
    \f[
    V = \pi G a_1 a_2 a_3 \int^\infty_0\frac{du}{\Delta}
    \left[\Psi(1) - \Psi(m^2(u))\right]
    \f]
    where
    \f[
    m^2(u) = \sum_{i=1}^3\frac{x_i^2}{a_i^2 + u}
    \f]
    and
    \f[
    \Psi(m^2) = \int_1^{m^2} dm^2\,\rho(m^2).
    \f]

    The external potential has the form (Cha ndra pg. 51, eqn. 89)
    \f[
    V = \pi G a_1 a_2 a_3 \int^\infty_\lambda\frac{du}{\Delta}
    \left[\Psi(1) - \Psi(m^2(u))\right].
    \f]
*/    
class EllipsoidForce 
{
public:

  enum BarType {powerlaw, ferrers, expon};

private:

  vector<double> a;
  double mass, param, rho0, rmin, rmax;
  BarType bartype;

  std::shared_ptr<LegeQuad> gq, gz, gt;

  int num, numr, N;

  vector<double> rr, mm, pp, uu;
  double lrmin, lrmax, ldr;

  double solve(vector<double> x, double m2);

  vector<int> ntab;
  vector<double> dX, table;
  double rtmin, rtmax;
  bool quadr, tmade, tlog;
  bool nindx(Eigen::Vector3d& x, Eigen::Vector3i& n);

  void write_cache();
  bool read_cache();

public:

  //! Number of points for the monopole/quadrupole computation
  static int n;

  //! Fraction of grid step for numerical derivatives
  static double dfrac;

  //! Interpolation grid cache name
  static string cache_file;

  //! Constructor
  EllipsoidForce(double a, double b, double c, double mass, 
		 double rmin, double rmax,
		 BarType bartype, bool quad, 
		 double param=0.0, int num=100, int numr=100);
  //! Destructor
  ~EllipsoidForce();

  //! Returns the total mass and moment of inertia vector
  void MassInertia(double& M, vector<double>& I);

  //! Returns the potential at position x
  double getPotl(vector<double> x);

  //! Returns the radial part of quadrupole potential
  double getPotl(double r);

  //! Returns the density at position x
  double getDens(vector<double> x);

  //! Compute the surface density along the major axis
  double getSurfDens(double r);

  //! Quadrupole bar density for radius r
  double RhoBar(double r);

  //! Quadrupole potential term for radius r
  double U22(double r);

  //! Get the quadrupole components at the radius r
  bool quadpot(double r, double& f, double& fr, double& M);

  //! Tabulate the potential and force
  void MakeTable(int n1, int n2, int n3);

  //! Evaluate the force (indices 0,1,2) potential (index 3) on the table
  void TableEval(vector<double> x, vector<double>& forcepot);
};

#endif // EllipsoidForce
