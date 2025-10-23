#ifndef _UserLogPot_H
#define _UserLogPot_H

#include "Particle.H"
#include "AxisymmetricBasis.H"
#include "ExternalCollection.H"

/** Logarithmic triaxial potential

    @param R is the core size
    @param b is the y-axis scale factor
    @param c is the z-axis scale factor
    @param v2 is the amplitude (e.g. the velocity on the axis at x=1)
*/
class UserLogPot : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double R, b, c, v2;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserLogPot(const YAML::Node& conf);

  //! Destructor
  ~UserLogPot();

};

#endif
