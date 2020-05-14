#ifndef _UserAddMass_H
#define _UserAddMass_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

#include <vector>
#include <array>

/** Add mass to component

    The particle positions are computing using either the spherical
    (Halo algorithm) or projected cylindrical (Disk algorithm) mass
    profile from the phase space.

    For the Disk algorithm, the particles are added to the mean disk
    midplane based on a real-time rotation curve estimate.  For the
    Halo algorithm, the velocity is chosen isotropically based on the
    mean velocity dispersion in spherical shells.

    @param compname is the component and must be specified
    @param tstart is the first time when particles are added
    @param dtime is the time interval between particle additions
    @param number is the number of particles to add per addition
    @param mass is the particle mass to add
    @param algorithm specifies the particle selection method
    @param rmin lower bin radius
    @param rmax upper bin radius
    @param numr is the number of radial bins
    @param logr true means use logarithmic bining
    @param seed for the pseudorandom number generator
**/

class UserAddMass : public ExternalForce
{
private:
  
  std::string comp_name;
  Component *c0;

  bool logr;
  unsigned number, numr;
  double mass, rmin, rmax, lrmin, lrmax, dr, tnext, dtime;
  long int seed;

  std::vector< std::vector<double> > mas_bins;
  std::vector< std::vector<std::array<double, 3>> > vl_bins, v2_bins;
  std::vector< std::vector<std::array<double, 3>> > L3_bins;

  void clear_bins();
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  void userinfo();

  enum class Algorithm { Disk, Halo } alg;

  static std::map<std::string, Algorithm> AlgMap;
  static std::map<Algorithm, std::string> AlgMapInv;

  boost::shared_ptr<ACG>     gen;
  boost::shared_ptr<Uniform> urand;
  boost::shared_ptr<Normal>  nrand;

public:

  //! Constructor
  UserAddMass(const YAML::Node& conf);

  //! Destructor
  ~UserAddMass();

};

#endif
