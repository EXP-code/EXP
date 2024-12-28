#ifndef _CBDisk_H
#define _CBDisk_H

#include <memory>
#include <map>

#include <PolarBasis.H>

/** Reference implemenation of the Clutton-Brock disk basis

    Parameters

    \param mmax is the maximum azimuthal order

    \param Lmax alias for mmax

    \param Mmax alias for mmax

    \param scale radius for coordinate scaling.  By default, the core radius of the Clutton-Brock density is 1.
*/
class CBDisk : public PolarBasis 
{

private:

  //@{
  //! 2D Clutton-Brock basis
  double phif(const int n, const int m, const double r);
  double dphi(const int n, const int m, const double r);

  double potl(const int n, const int m, const double r);
  double dpot(const int n, const int m, const double r);
  double dens(const int n, const int m, const double r);

  double norm(const int n, const int m);

  void potl(const int m, const double r, Eigen::VectorXd& a);
  void dens(const int m, const double r, Eigen::VectorXd& a);
  //@}
  
  //! Required maximum table radious (no limit here)
  virtual double getRtable() { return std::numeric_limits<double>::max(); }

  //! Sanity check
  void orthoCheck();

  void initialize(void);

  void get_dpotl(double r, double z,
		 Eigen::MatrixXd& p, Eigen::MatrixXd& dpr, Eigen::MatrixXd& dpz, int tid);

  void get_potl(double r, double z, Eigen::MatrixXd& p, int tid);

  //! Background instance
  using Disk2d = std::shared_ptr<EmpCyl2d::ModelCyl>;
  Disk2d disk;

  //! Set background from YAML
  void setBackground();

  //! Background evaluation
  virtual std::tuple<double, double, double>
  background(double r, double z)
  { return {disk->pot(r), -disk->dpot(r), 0.0}; }
  
  void get_dens(double r, double z, Eigen::MatrixXd& p, int tid);

  void get_potl_dens(double r, double z,
		     Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid);

  //@{
  //! Parameters
  double scale;
  int    mmax;
  string model;
  //@}

  //! Potential, force, and density scale factors
  double fac1, fac2;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

protected:

  //! Clutton-Brock potential
  void get_pot(Eigen::MatrixXd& Vc, Eigen::MatrixXd& Vs, double r, double z)
  {
    get_potl(r, z, Vc, 0);
    Vs = Vc;
  }

public:
				// Global parameters
  /** Constructor
      @param c0 is the instantiating caller (Component)
      @param conf passes any parameters to basis instance
      @param scale is the radius for coordinate scaling
  */
  CBDisk(Component* c0, const YAML::Node& conf, MixtureBasis* m=0);

  //! Destructor
  virtual ~CBDisk();
};

#endif


