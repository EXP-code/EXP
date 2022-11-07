#ifndef _UserSlabHalo_H
#define _UserSlabHalo_H

#include <Particle.H>
#include <ExternalCollection.H>

/** Sech^2 external 1-dimensional halo model

    @param ctrname component will be used to define the center
    @param h0 is the slab scale height
    @param rho0 is the slab central density
    @param v0 is the one-dimensional isothermal velocity dispersion
    @param z0 is the midplane position
*/
class UserSlabHalo : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  string ctr_name;
  double h0, rho0, sig20, U0, v0, v20, z0;
  Component *c0;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserSlabHalo(const YAML::Node& conf);

  //! Destructor
  ~UserSlabHalo();

};

#endif
