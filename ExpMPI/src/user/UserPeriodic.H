#ifndef _UserPeriodic_H
#define _UserPeriodic_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Enforce periodic boundary conditions
    @param ctrname is the name of the component whose particles will be subject to the boundary condition
    @param sx is the cube size in the x-direction
    @param sy is the cube size in the y-direction
    @param sz is the cube size in the z-direction
    @param cx is the center offset in the x-direction
    @param cy is the center offset in the y-direction
    @param cz is the center offset in the z-direction
    @param btype is a string of length that defines the type of boundary condition for each dimension ('r' for reflection bc and 'p' for periodic bc)
*/
class UserPeriodic : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  vector<double> offset, L;
  string comp_name, bc;

  Component *c0;

  void userinfo();

public:

  //! Constructor
  UserPeriodic(string &line);

  //! Destructor
  ~UserPeriodic();

};

#endif
