#ifndef _UserDiffuse_H
#define _UserDiffuse_H

#include <random>

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <global.H>

/** Random walk of each orbit on energy surface

    @param rate of diffusion
    @param seed for the pseudorandom number generator
    @param nfreq steps between diffusion coefficient recomputation
    @param rmin is the minimum radiuis of coefficient grid
    @param rmax is the maximum radiuis of coefficient grid
    @param pmass is the mass of virtual two-body particle
    @param logL for diffusion coefficients
    @param clip the Gaussian variates at this number of "Sigma"
    @param numr is the number of radial entries in diff. coefficient table
    @param numv is the number of velocity  entries in diff. coefficient table
    @param nume is the number of energy integration knots for diff. coef.
    @param logr is set to 1 for radial logarithmic spacing in diff. coef. table
    @param use_file is set to 1 to read background model from "modfile"
    @param check_ev is set to 1 to check coordinate rotations (debug only; don't use this for production)
    @param modfile contains the background model
    @param diverge set to extrapolate model for a power law cusp
    @param diverge_rfac exponent for the power law
*/
class UserDiffuse : public ExternalForce
{

  typedef std::vector<double> dvector;

private:
  
  string name;
  long int seed;
  int nfreq;
  double rmin, rmax;
  double pmass, logL, clip;
  int numr;
  int numv;
  int nume;
  bool logr;
  bool use_file;
  bool check_ev;
  string modfile;
  int diverge;
  double diverge_rfac;

  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();
  void userinfo();

  SphericalModelTable *model;

  void compute_model();
  void compute_diffuse();
  void get_coefs(double r, double v,
		 double& dvpara1, double& dvpara2, double& dvperp2);


  double delR, delV, Gamma, Rmin, Rmax, Emin, Emax;
  vector<double> R, Vmax;
  vector<dvector> dVpara1, dVpara2, dVperp2;

  LegeQuad *jq;

  std::uniform_real_distribution<> urand;
  std::normal_distribution<> nrand;

  string check_file;
  std::vector< std::vector<double> > ev_mean_th, ev_disp_th;
  std::vector<int> ev_numb_th;

public:

  //! Constructor
  UserDiffuse(const YAML::Node& conf);

  //! Destructor
  ~UserDiffuse();

};

#endif
