#ifndef _Basis_H
#define _Basis_H

#include <PotAccel.H>
#include <Eigen/Eigen>

//! Defines a basis-based potential and acceleration class
class Basis : public PotAccel
{
protected:

public:

  //! Constructor
  Basis(Component* c0, const YAML::Node& conf);

  //! All derived classes must provide a method to compute the force
  virtual void determine_acceleration_and_potential(void) = 0;

  /** Supply density, potential and derivatives at a point in
      Cartesian coordinates */
  virtual void determine_fields_at_point
  (double x, double y, double z,
   double *tdens0, double *tpotl0, double *tdens, double *tpotl, 
   double *tpotx, double *tpoty, double *tpotz) = 0;

  /** Supply density, potential and derivatives at a point in
      spherical coordinates */
  virtual void determine_fields_at_point_sph
  (double r, double theta, double phi,
   double *tdens0, double *dpotl0, double *tdens, double *tpotl, 
   double *tpotr, double *tpott, double *tpotp) = 0;


  /** Supply density, potential and derivatives at a point in
      cylindrical coordinates */
  virtual void determine_fields_at_point_cyl
  (double r, double z, double phi,
   double *tdens0, double *dpotl0, double *tdens, double *tpotl,
   double *tpotr, double *tpotz, double *tpotp) = 0;


  /** @name Utility functions */
  // @{

  //! Compute Associated Legendre Polynomials, return MatrixXd type
  void legendre_R(int lmax, double x, Eigen::MatrixXd& p);

  /** Compute Associated Legendre Polynomials and derivitives, 
  return Matrix type */
  void dlegendre_R(int lmax, double x, Eigen::MatrixXd &p, Eigen::MatrixXd &dp);

  //! Compute vectors of sines and cosines by recursion
  void sinecosine_R(int mmax, double phi, Eigen::VectorXd& c, Eigen::VectorXd& s);

  // @}

};

#endif // _Basis_H


