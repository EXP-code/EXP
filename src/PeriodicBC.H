#ifndef _PeriodicBC_H
#define _PeriodicBC_H

#include <random>

#include "Particle.H"
#include "AxisymmetricBasis.H"
#include "ExternalCollection.H"

/** Enforce periodic boundary conditions in a rectangular prism (box)
    @param ctrname is the name of the component whose particles will be subject to the boundary condition
    @param sx is the edge size in the x-direction
    @param sy is the edge size in the y-direction
    @param sz is the edge size in the z-direction
    @param cx is the center offset in the x-direction
    @param cy is the center offset in the y-direction
    @param cz is the center offset in the z-direction
    @param btype is a string of length that defines the type of boundary condition for each dimension ('r' for reflection bc, 'p' for periodic bc, 'v' is vacuum)
*/
class PeriodicBC : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  std::vector<double> offset, L;
  std::string comp_name, bc;

  // Needed for converting thermal energy using Boltzmann constant
  double vunit;
  double munit;
  double temp;

  std::vector<double> atomic_weights;

  Component *c0;

  void userinfo();

  //! Pseudorandom generation
  std::normal_distribution<> norm;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  PeriodicBC(const YAML::Node& conf);

  //! Destructor
  ~PeriodicBC();

};

#endif
