#ifndef _UserAtmos_H
#define _UserAtmos_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Linear gravitational field

    @param gx component of the acceleration in the x-direction
    @param gy component of the acceleration in the y-direction
    @param gz component of the acceleration in the z-direction
    @param compname is the name of the component accelerated by the field
*/
class UserAtmos : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  vector<double> g;
  string compname;

  Component *c0;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserAtmos(const YAML::Node& conf);

  //! Destructor
  ~UserAtmos();

};

#endif
