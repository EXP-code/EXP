#ifndef _Hernquist_H
#define _Hernquist_H

#include <SphericalBasis.H>

class MixtureBasis;

//! Biorthonormal pairs for the Hernquist model
class Hernquist : public SphericalBasis
{

private:

  void get_dpotl(int lmax, int nmax, double r,
		 Eigen::MatrixXd& p, Eigen::MatrixXd& dp, int tid);
  void get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);
  void get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);
  void get_potl_dens(int lmax, int nmax, double r,
		     Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid);

  double norm(int, int);
  double knl(int, int);

  double xi_to_r  (double x);
  double d_r_to_xi(double r);
  double r_to_xi  (double r);

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  Hernquist(Component* c0, const YAML::Node& conf, MixtureBasis* m=0);

};

#endif
