#ifndef _UserMassEvo_H
#define _UserMassEvo_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Evolve the per particle masses for given component using a cubic polynomial

    The value t0 is assumed to be the initial time of the simulation.
    Then, the mass per particle is evolved according to:

    m *= a0 + a1*(t - t0) + a2*(t - t0)**2 + a3*(t - t0)**3

    @param comp_name is the component and must be specified
    @param t0 is the zero time (0 by default)
    @param a0 is the constant term (1 by default)
    @param a1 is the linear term (0 by default)
    @param a2 is the quadratic term (0 by default)
    @param a3 is the cubic term (0 by default)

    For example, to linearly increase the mass of a particle by a
    factor of 10 over 5 time units use:
    t0 = 0
    a0 = 1
    a1 = 2.5
    a2 = a3 = 0
*/
class UserMassEvo : public ExternalForce
{
private:
  
  std::string comp_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  //! Parameters
  double t0, a0, a1, a2, a3;

  //! Current value of mass-evolution function
  double curr;

  //! Last value of mass-evolution function
  double last;

  void userinfo();

  double fct(double t)
  {
    double dt = t - t0;
    return a0 + dt*(a1 + dt*(a2 + a3*dt));
  }

public:

  //! Constructor
  UserMassEvo(const YAML::Node& conf);

  //! Destructor
  ~UserMassEvo();

};

#endif
