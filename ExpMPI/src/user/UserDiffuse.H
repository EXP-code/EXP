#ifndef _UserDiffuse_H
#define _UserDiffuse_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Random walk of each orbit on energy surface

    @param rate of diffusion
    @param seed for the pseudorandom number generator
*/
class UserDiffuse : public ExternalForce
{
private:
  
  string name;
  long int seed;
  double rate; 

  Component *c0;

  void determine_acceleration_and_potential(void);
  void determine_acceleration_and_potential(vector<Particle>*);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();
  void userinfo();

  ACG *gen;
  Normal *nrand;

public:

  //! Constructor
  UserDiffuse(string &line);

  //! Destructor
  ~UserDiffuse();

};

#endif
