#ifndef _UserHalo_H
#define _UserHalo_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Halo model based on SphericalModelTable

    @param model_file is the filename
    @param diverge set to 1 uses a power law inner profile
    @param diverge_rfac is the value of power law exponent
*/
class UserHalo : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void determine_acceleration_and_potential(vector<Particle>*);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  string model_file;
  int diverge;
  double diverge_rfac;

  SphericalModelTable *model;

  void userinfo();

public:

  //! Constructor
  UserHalo(string &line);

  //! Destructor
  ~UserHalo();

};

#endif
