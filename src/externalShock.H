#ifndef _externalShock_H
#define _externalShock_H

#include <orbit.h>
#include <massmodel.h>

//! Compute shock perturbation
class externalShock : public ExternalForce
{

private:

  double E, K, PER, AMPL;
  string INFILE;

  SphericalModelTable *model;
  SphericalOrbit *t;

  double get_tidal_shock(double T);

  void initialize();

  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! Constructor
  externalShock(const YAML::Node& conf);

  //! Destructor
  ~externalShock();

};


#endif

