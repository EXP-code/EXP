#ifndef _UserMNdisk_H
#define _UserMNdisk_H

#include "Particle.H"
#include "AxisymmetricBasis.H"
#include "ExternalCollection.H"

/** Potential and force of the Miyamoto-Nagai disk.

    @param a is the disk scale length parameter
    @param b is the disk scale height parameter
    @param mass is the total mass of the disk
    @param Ton is the time at the center of the <code>erf</code> turn on
    @param Toff is the time at the center of the <code>erf</code> turn off
    @param DeltaT is the spread of the turn on
    @param ctrname defines the component that defines the center (otherwise the center is fixed to the origin)
*/
class UserMNdisk : public ExternalForce
{
private:
  
  string ctr_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double a, b, mass, Ton, Toff, DeltaT;

  string name;

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserMNdisk(const YAML::Node& conf);

  //! Destructor
  ~UserMNdisk();

};

#endif
