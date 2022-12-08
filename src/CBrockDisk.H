#ifndef _CBrockDisk_H
#define _CBrockDisk_H

#include <memory>
#include <vector>
#include <string>
#include <set>

#include <AxisymmetricBasis.H>
#include <CoefContainer.H>
#include <TwoDCoefs.H>

class MixtureBasis;

//! This routine computes the potential, acceleration and density
//! using the Clutton-Brock flat disk expansion
class CBrockDisk : public AxisymmetricBasis
{

private:

  MixtureBasis *mix;
  pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;
  double *cylmass0;

  std::shared_ptr<Eigen::MatrixXd> expcoef, expcoef1, expcoefP;
  std::vector<Eigen::MatrixXd> expcoef0, cc, cc1;
  Eigen::MatrixXd normM, dend, work;
  std::vector<Eigen::MatrixXd> potd, dpot;
  std::vector<Eigen::VectorXd> u, du, cosm, sinm;

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
  void determine_coefficients_playback();
  void determine_coefficients_particles();
  void determine_acceleration_and_potential();

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

				// Parameters
  double rmax, scale;
  int Lmax, nmax;
  bool self_consistent;

  //! Coefficient magic number
  const unsigned int cmagic = 0xc0a57a4;

  //! Playback object
  std::shared_ptr<TwoDCoefs> playback;

  /** Master node ships coefficients to hosts.  True (default) implies
      that only the master node caches the coefficients for playback
      to save core memory.  This is set in the config input using the
      'coefMaster: bool' parameter. Once I am sure that there are no
      algorithmic issues, I will remove this as an option.
  */
  bool coefMaster;

  //! Last playback coefficient evaluation time
  double lastPlayTime;

  //! Coefficient container instance for writing HDF5
  CoefClasses::CylCoefs cylCoefs;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  //! \param c0 is the instantiating caller (Component)
  //! \param line is the parameter string
  //! \param m is the MixtureBasis for a multicenter expansion
  CBrockDisk(Component* c0, const YAML::Node& conf, MixtureBasis* m=0);

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

  //! Save coefficients to HDF5 file
  void dump_coefs_h5(const std::string& file);
};


#endif
