#ifndef _UserAddMass_H
#define _UserAddMass_H

#include <random>

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <global.H>

#include <vector>
#include <array>
#include <list>

/** Add mass to component

    The particle positions are computing using either the spherical
    (Halo algorithm) or projected cylindrical (Disk algorithm) mass
    profile from the phase space.

    For the Disk algorithm, the particles are added to the mean disk
    midplane based on a real-time rotation curve estimate.  For the
    Halo algorithm, the velocity is chosen isotropically based on the
    mean velocity dispersion in spherical shells.  

    If the accel parameter is true, we assign the initial acceleration
    for created particles from the basis force fields.

    This routine should occur in the [External] sequence before any
    true force evalaution to include the external force.

    @param compname is the component and must be specified
    @param complist is the list of components for force evaluation, if specified
    @param tstart is the first time when particles are added
    @param dt is the time interval between particle additions
    @param number is the number of particles to add per addition
    @param mass is the particle mass to add
    @param algorithm specifies the particle selection method
    @param rmin lower bin radius
    @param rmax upper bin radius
    @param numr is the number of radial bins
    @param logr true means use logarithmic bining
    @param seed for the pseudorandom number generator
    @param debug is the time interval between profile outputs (<=0 for off)
    @param interp set to true for interpolation in angular momentum bins
    @param planar set to true enforces particles added in the x-y plane for testing (default: false)
    @param accel set to true assigns the acceration by evaluating the basis force (default: true)
    @param maxtry is the maximum number of attempts to select a particle from the mass array
**/
class UserAddMass : public ExternalForce
{
private:

  std::vector<std::string> comp_list;
  std::string comp_name, fvdisp;
  std::vector<Component*> comp_vec;
  list<Component*> comp_acc;
  Component *c0;

  bool logr, interp, cforce, planar, accel;
  unsigned number, numr;
  double mass, rmin, rmax, lrmin, lrmax, dr, tnext, dt, debug, dnext;
  double tstart, scale, vdisp;
  double Rd, drprof, Zd;
  long int seed, numRprof, numzprof, numvdisp;
  int maxtry;
  
  std::vector< std::vector<double> > mas_bins_current;
  std::vector< double > Mprof, Rprof, vdispR, vdispV;


  const double ttol = 1.0e-8;


  using vec3 = std::array<double, 3>;
  std::vector< std::vector<vec3> > vl_bins, v2_bins;
  std::vector< std::vector<vec3> > L3_bins;
  std::vector< std::vector<double> > z_mean_bins, vz_mean_bins;


  void clear_bins();
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  void userinfo();

  enum class Algorithm { Disk, Halo } alg;

  static std::map<std::string, Algorithm> AlgMap;
  static std::map<Algorithm, std::string> AlgMapInv;

  std::uniform_real_distribution<> urand;
  std::normal_distribution<> nrand;

public:

  //! Constructor
  UserAddMass(const YAML::Node& conf);

  //! Destructor
  ~UserAddMass();

};

#endif
