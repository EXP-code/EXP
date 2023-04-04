#ifndef _UserAgnNoise_H
#define _UserAgnNoise_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Evolve the per particle masses after after AGN/star burst event

    @param compname is the component and must be specified
    @param tau is the half-life between events (0.1 by default)
    @param R0 is the radius of the event (0.003 by default)
    @param eps is the amplitude of the event (0.1 by default)
*/
class UserAgnNoise : public ExternalForce
{
private:
  
  std::string comp_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  //@{
  //! Parameters
  double tau, R0, eps;
  int    loc;
  //@}

  //! Current value of mass-evolution function
  double tev;

  void userinfo();

  //! For selecting events
  std::uniform_real_distribution<double> number_01;

public:

  //! Constructor
  UserAgnNoise(const YAML::Node& conf);

  //! Destructor
  ~UserAgnNoise();

};

#endif
