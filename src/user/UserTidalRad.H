#ifndef _UserTidalRad_H
#define _UserTidalRad_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Compute tidal radius of a component using a total energy criterion

    The radius is estimated by computing the mean radius of particles
    in a fixed percentile of the distribution between the maximally
    and minimally bound particle ranked by binding energy.  The
    average radius is then multiplied by a factor of the original
    value from the measured percentile to the tidal or scaling radius
    for the estimate.  In tests, a percentile of 0.9 seems to be a
    stable choice.  Clearly, choosing the minimally bound particle as
    the estimator is way too noisy.

    @param comp_name is the component and must be specified
    @param rtrunc is the initial truncation radius
    @param rfactor multiplies the computed truncation radius for setting 
    @param rtorig is the target percentile radius at simulation start
    @param pctile is the target percentile for radial average
    @param pcnbin is the number of particles around the percentile bin
    @param boxcar is the number of time steps for radial average
    @param dtTrunc is the time between truncation radius updates (0 means never)
    @param dtScale is the time between scale radius updates (0 means never)
    @param debug set to true prints mean percentile radii for examination
    @param diag is the stride for diagnostic output (0 means never)
*/
class UserTidalRad : public ExternalForce
{
private:
  
  std::string comp_name;
  double rtrunc, rfactor, pctile, dtTrunc, dtScale, rtorig;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  std::string filename;
  std::vector<double> cov, mas, erg, rad;
  std::vector<double> radsol;
  double rt_cur, tnextT, tnextS;
  int pass, pcnbin, diag;
  unsigned boxcar;
  bool debug;
  
  void userinfo();

  double radavg();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserTidalRad(const YAML::Node& conf);

  //! Destructor
  ~UserTidalRad();

};

#endif
