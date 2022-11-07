#ifndef _CBrock_H
#define _CBrock_H

#include <SphericalBasis.H>

//! Biorthonormal pairs based on Clutton-Brock's ultraspherical series
class CBrock : public SphericalBasis
{

private:

  Eigen::VectorXd uu, duu;

  void get_dpotl(int lmax, int nmax, double r,
		 Eigen::MatrixXd& p, Eigen::MatrixXd& dp, int tid);
  void get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);
  void get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);
  void get_potl_dens(int lmax, int nmax, double r,
		     Eigen::MatrixXd& p, Eigen::MatrixXd& d,int tid);

  void initialize(void);

  double knl(int, int);
  double norm(int, int);

  double r_to_rh(double r);
  double rh_to_r(double rh);
  double d_r_to_rh(double);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  /*! \param c0 is the instantiating component */
  /*! \param line is the parameter string */
  /*! \param m is the MixtureBasis for a multicenter expansion */
  CBrock(Component* c0, const YAML::Node& conf, MixtureBasis* m=0);

};

#endif
