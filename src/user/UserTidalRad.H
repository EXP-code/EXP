#ifndef _UserTidalRad_H
#define _UserTidalRad_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Compute tidal radius of a component using a total energy criterion

    The radius is estimated by computing the mean radius of particles
    in a fixed percentile of the distribution between the maximally
    and minimally bound particle ranked by binding energy.  The
    average radius is then multiplied by a factor (of order 2) for the
    estimate.  In tests, a percentile of 0.9 with a factor 2 seems to
    be a stable choice.  Clearly, choosing the minimally bound
    particle as the estimator is way too noisy.

    @param comp_name is the component and must be specified
    @param rtrunc is the initial truncation radius
    @param rfactor multiplies the computed truncation radius for setting 
    @param dryrun set to true does not apply basis scaling
    @param pctile is the target percentile for radial average
    @param pcnbin is the number of particles around the percentile bin
    @param boxcar is the number of time steps for radial average
    @param update is the number of time steps between scale updates
    @param debug set to true prints mean percentile radii for examination
*/
class UserTidalRad : public ExternalForce
{
private:
  
  std::string comp_name;
  double rtrunc, rfactor, pctile;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  std::string filename;
  std::vector<double> cov, mas, erg, rad;
  std::vector<double> radsol;
  double rt_cur, rt_cur0;
  bool firsttime, dryrun, debug;
  int pass, pcnbin, update;
  unsigned boxcar;
  
  void userinfo();

  double radavg();

public:

  //! Constructor
  UserTidalRad(const YAML::Node& conf);

  //! Destructor
  ~UserTidalRad();

};

#endif
