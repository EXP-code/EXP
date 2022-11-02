#ifndef _UserShear_H
#define _UserShear_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Shearing sheet equations of motion based on gravitational
    potential for a flat rotation curve.

    @param radius is the value of R_0 in radial units
    @param velocity times sqrt(2) is the circular velocity in velocity units
*/
class UserShear : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  string ctr_name;
  double r0, s0, omega, kappa, A, B, xoff;
  Component *c0;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserShear(const YAML::Node& conf);

  //! Destructor
  ~UserShear();

};

#endif
