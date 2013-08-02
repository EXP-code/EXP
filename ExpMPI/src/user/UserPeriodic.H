#ifndef _UserPeriodic_H
#define _UserPeriodic_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Enforce periodic boundary conditions in a rectangular prism (box)
    @param ctrname is the name of the component whose particles will be subject to the boundary condition
    @param sx is the cube size in the x-direction
    @param sy is the cube size in the y-direction
    @param sz is the cube size in the z-direction
    @param cx is the center offset in the x-direction
    @param cy is the center offset in the y-direction
    @param cz is the center offset in the z-direction
    @param btype is a string of length that defines the type of boundary condition for each dimension ('r' for reflection bc and 'p' for periodic bc)
    @param nbin is the number of bins in the x-direction (0 signals no trace)
    @param dT is the time step for binned density/temp trace
    @param tcol is the double attribute index for the temperature data (ignored if < 0)
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

  vector< vector<unsigned> > cbinT;
  vector< vector<double>   > mbinT, tbinT;
  vector<unsigned>           cbin;
  vector<double>             mbin, tbin;
  bool trace;

  Component *c0;

  void userinfo();

public:

  //! Constructor
  UserPeriodic(string &line);

  //! Destructor
  ~UserPeriodic();

};

#endif
