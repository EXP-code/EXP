#ifndef _UserPeriodic_H
#define _UserPeriodic_H

#include <random>

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Enforce periodic boundary conditions in a rectangular prism (box)
    @param ctrname is the name of the component whose particles will be subject to the boundary condition
    @param sx is the edge size in the x-direction
    @param sy is the edge size in the y-direction
    @param sz is the edge size in the z-direction
    @param cx is the center offset in the x-direction
    @param cy is the center offset in the y-direction
    @param cz is the center offset in the z-direction
    @param btype is a string of length that defines the type of boundary condition for each dimension ('r' for reflection bc, 'p' for periodic bc, 'v' is vacuum, and 't' is thermal)
    @param nbin is the number of bins in the x-direction (0 signals no trace)
    @param dT is the time step for binned density/temp trace
    @param tcol is the double attribute index for the temperature data (ignored if < 0)
    @param temp is the temperature of the wall for 'thermal' boundary conditions (ignored for r, p, and v boundaries)
*/
class UserPeriodic : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  void write_trace();

  vector<double> offset, L;
  string comp_name, bc;
  double dT, Tnext, dX;
  int nbin, tcol;
  bool thermal;

  // Needed for converting thermal energy using Boltzmann constant
  double vunit;
  double munit;
  double temp;

  std::vector<double> atomic_weights;

  vector< vector<unsigned> > cbinT;
  vector< vector<double>   > mbinT, tbinT;
  vector<unsigned>           cbin;
  vector<double>             mbin, tbin;
  bool trace;

  Component *c0;

  void userinfo();

  //! Pseudorandom generation
  std::normal_distribution<> norm;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserPeriodic(const YAML::Node& conf);

  //! Destructor
  ~UserPeriodic();

};

#endif
