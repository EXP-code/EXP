#ifndef _UserReflect_H
#define _UserReflect_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

#include <ACG.h>
#include <Uniform.h>

/** Enforce reflecting boundary conditions for sphere
    @param compname is the name of the component whose particles will be subject to the boundary condition
    @param radius is the radius of the reflecting boundary
*/
class UserReflect : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  void write_trace();

  string comp_name;
  double radius;
  bool debug;
  unsigned *wrong_dir, *too_big;
  vector<ACG*> gen;
  vector<Uniform*> unit;

  Component *c0;

  void userinfo();

public:

  //! Constructor
  UserReflect(const YAML::Node& conf);

  //! Destructor
  ~UserReflect();

};

#endif
