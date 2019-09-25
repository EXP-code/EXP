#ifndef _UserTidalRad_H
#define _UserTidalRad_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Compute tidal radius of a component using a total energy criterion

    @param comp_name is the component and must be specified
    @param rtrunc is the initial truncation radius
*/
class UserTidalRad : public ExternalForce
{
private:
  
  std::string comp_name;
  double rtrunc, rfactor;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  std::string filename;
  std::vector<double> rad, erg, cov, mas;
  double rt_cur, rt_cur0;
  bool firsttime;
  int pass;
  
  void userinfo();

public:

  //! Constructor
  UserTidalRad(const YAML::Node& conf);

  //! Destructor
  ~UserTidalRad();

};

#endif
