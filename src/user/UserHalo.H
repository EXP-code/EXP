#ifndef _UserHalo_H
#define _UserHalo_H

#include "Particle.H"
#include "AxisymmetricBasis.H"
#include "ExternalCollection.H"

/** Halo model based on SphericalModelTable

    @param model_file is the filename
    @param comp_name is the component to apply this (blank for all, default="")
    @param diverge set to 1 uses a power law inner profile
    @param diverge_rfac is the value of power law exponent
    @param q0 is the x-axis scaling (for simulating a triaxial halo, default: 1)
    @param q1 is the y-axis scaling (for simulating a triaxial halo, default: 1)
    @param q2 is the z-axis scaling (for simulating a triaxial halo, default: 1)
*/
class UserHalo : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  std::string model_file, comp_name;
  double q1, q2, q3;
  double diverge_rfac;
  int diverge;

  SphericalModelTable *model;
  Component *c0;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserHalo(const YAML::Node& conf);

  //! Destructor
  ~UserHalo();

};

#endif
