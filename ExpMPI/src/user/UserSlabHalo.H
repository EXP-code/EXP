#ifndef _UserSlabHalo_H
#define _UserSlabHalo_H

#include <Particle.H>
#include <ExternalCollection.H>

/** Sech^2 external 1-dimensional halo model

    @param model_file is the filename
    @param diverge set to 1 uses a power law inner profile
    @param diverge_rfac is the value of power law exponent
*/
class UserSlabHalo : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  string ctr_name;
  double h0, rho0, sig20, U0, v0, v20;
  Component *c0;

  void userinfo();

public:

  //! Constructor
  UserSlabHalo(string &line);

  //! Destructor
  ~UserSlabHalo();

};

#endif
