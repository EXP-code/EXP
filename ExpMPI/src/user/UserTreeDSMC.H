#ifndef _UserTreeDSMC_H
#define _UserTreeDSMC_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

#include <Timer.h>

#include "CollideLTE.H"

/** DSMC particle routine
    @param Lunit is the physical scale in cm per system units
    @param Tunit is the physical time scale in years per system unit
    @param Munit is the physical mass in units of mp per system units
    @param BoxSize half the size of a side in system units
    @param ncell is the desired number of particles per cell
    @param cnum non-zero adjusts cross section for cnum collisions per cel
    @param diamfac is a multiplicative factor for the Bohr radius (converted to system units
    @param taufrac is a multiplicative fraction of the crossing time for time step selection
*/
class UserTreeDSMC : public ExternalForce
{
private:
  
  int ncell;
  int cnum;
  int nsteps;
  double diamfac;
  double taufrac;
  double boxsize;


  CollideLTE *collide;

  string comp_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

				// The tree
  pHOT tree;
  double volume;

				// Timers
  Timer driftTime, bcTime, partnTime, treeTime;
  Timer collideTime, densityTime, gatherTime;

  void userinfo();

#ifdef TIMER
  Timer *timer;
#endif

public:

				// System to physical scaling factors
  static double Lunit;
  static double Tunit;
  static double Vunit;
  static double Munit;
  static double Eunit;

  //! Constructor
  UserTreeDSMC(string &line);

  //! Destructor
  ~UserTreeDSMC();
};

#endif
