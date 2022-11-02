#ifndef _UserBarTrigger_H
#define _UserBarTrigger_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Mulitpole external force computed from satellite interlopers on
    linear trajectories to be used explore bar triggering

    @param ctrname is the component to center upon
    @param smass is the mass of satellite in system units
    @param impact is the impact parameter in system units
    @param theta is the orbital plane inclination in degrees
    @param svel is the velocity of the satellite in system units
    @param stime is the offset of the satellite where T=0 is the point
    of closest approach.  That is, stime=-3 implies that the satellite
    is svel*stime away from the point of closest approach when the
    system time is 0
    @param lmax is the maximum number of harmonics to keep
*/
class UserBarTrigger : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  string ctr_name;
  double impact, theta, smass, svel, stime;
  int lmax;

  Component *c0;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserBarTrigger(const YAML::Node& conf);

  //! Destructor
  ~UserBarTrigger();

};

#endif
