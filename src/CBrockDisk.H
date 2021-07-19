#ifndef _CBrockDisk_H
#define _CBrockDisk_H

#include <vector>
#include <AxisymmetricBasis.H>

class MixtureBasis;

//! This routine computes the potential, acceleration and density using the Clutton-Brock flat disk expansion
class CBrockDisk : public AxisymmetricBasis
{

private:

  MixtureBasis *mix;
  pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;
  double *cylmass0;

  Eigen::MatrixXd expcoef, expcoef1;
  vector<Eigen::MatrixXd> expcoef0, cc, cc1;
  Eigen::MatrixXd normM, dend, work;
  vector<Eigen::MatrixXd> potd, dpot;
  vector<Eigen::VectorXd> u, du, cosm, sinm;

  int use1, use0;

  void get_dpotl(int lmax, int nmax, double r, Eigen::MatrixXd& p, Eigen::MatrixXd& dp);
  void get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p);
  void get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p);
  void get_potl_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, Eigen::MatrixXd& d);
  double norm(int,int);

  void get_pot_coefs(int, const Eigen::VectorXd&, double &, double &);
  void get_dens_coefs(int, const Eigen::VectorXd&, double &);

  void get_pot_coefs_safe(int, const Eigen::VectorXd&, double &, double &, Eigen::MatrixXd& p, Eigen::MatrixXd& d);

  void initialize(void);

				// These can be hidden from the interface
  void determine_coefficients();
  void determine_acceleration_and_potential();

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

				// Parameters
  double rmax, scale;
  int Lmax, nmax;
  bool self_consistent;

public:

  //! Constructor
  //! \param line is the parameter string
  //! \param m is the MixtureBasis for a multicenter expansion
  CBrockDisk(const YAML::Node& conf, MixtureBasis* m=0);

  //! Destructor
  virtual ~CBrockDisk(void);

  //! Force evaluation called clients
  void get_acceleration_and_potential(Component*);


  /**
     Supply density, potential and derivatives at arbitrary point
     in polar coordinates.  Currently this is used to implement the
     required <code>determine_fields_at_point</code> members by
     ignoring the <code>z</code> coordinate.
  */
  void 
  determine_fields_at_point(double x, double y, double z,
			    double *tdens0, double *tpotl0, 
			    double *tdens, double *tpotl, 
			    double *tpotX, double *tpotY, 
			    double *tpotZ);

  void 
  determine_fields_at_point_polar(double r, double phi,
				  double *tdens0, double *tpotl0, 
				  double *tdens, double *tpotl, 
				  double *tpotr, double *tpotp);
  void 
  determine_fields_at_point_sph(double r, double theta, double phi,
				double *tdens0, double *tpotl0, 
				double *tdens, double *tpotl, 
				double *tpotr, double *tpott, 
				double *tpotp);

  void 
  determine_fields_at_point_cyl(double r, double z, double phi,
				double *tdens0, double *tpotl0, 
				double *tdens, double *tpotl, 
				double *tpotr, double *tpotz, 
				double *tpotp);

  //! Save coefficients to file (need type marker to id dump, component id?)
  void dump_coefs(ostream& out);

};


#endif
