#ifndef _UserDiffRot_H
#define _UserDiffRot_H

#include <Particle.H>
#include <Component.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Random walk of each orbit on energy surface

    @param rate of diffusion
    @param seed for the pseudorandom number generator
*/
class UserDiffRot : public ExternalForce
{
private:
  
  string name, avoid;
  long int seed;
  int maxpm;
  double rate; 
  bool first;

  Component *c0, *c1;
  vector<double> pos;
  int ipm, indx;

  void determine_acceleration_and_potential(void);
  void determine_acceleration_and_potential(vector<Particle>*);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();
  void userinfo();

  double get_dtime(Particle& p);

  ACG *gen;
  Uniform *uniform;

public:

  //! Constructor
  UserDiffRot(string &line);

  //! Destructor
  ~UserDiffRot();

};

#endif
