#ifndef _UserAgnNoise_H
#define _UserAgnNoise_H

#include <expand.H>
#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Evolve the per particle masses after after AGN/star burst event

    @param compname is the component and must be specified
    @param tau1 is the half-life between events (0.1 by default)
    @param tau2 is the recovery time (0.05 by default)
    @param R0 is the radius of the event (0.003 by default)
    @param eps is the amplitude of the event (0.1 by default)
*/
class UserAgnNoise : public ExternalForce
{
private:
  
  std::string comp_name;
  Component *c0;

#if HAVE_LIBCUDA==1
  //! Cuda implementation
  void determine_acceleration_and_potential_cuda();
#endif
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  //@{
  //! Parameters
  double tau1, tau2, R0, eps;
  int    loc;
  bool   info;
  //@}

  //! Current value of next AGN/starburst event
  double tev;

  //! For selecting events
  std::uniform_real_distribution<double> number_01;

  //! Log parameter info
  void userinfo();

public:

  //! Constructor
  UserAgnNoise(const YAML::Node& conf);

  //! Destructor
  ~UserAgnNoise();

};

#endif
