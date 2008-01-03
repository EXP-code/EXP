#ifndef _UserTreeDSMC_H
#define _UserTreeDSMC_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

#include <Timer.h>

#include "CollideLTE.H"

//! Helper class to pass info to threaded member
struct dsmc_pass_arguments
{
  //! Pointer to iterator
  pHOT *tree;
  //! Thread counter id
  int id;
};
  
class UserTreeDSMC;

struct thrd_pass_dsmc
{
  UserTreeDSMC *p;
  //! Arguments
  dsmc_pass_arguments arg;
};


/** DSMC particle routine
    @param Lunit is the physical scale in cm per system units
    @param Tunit is the physical time scale in years per system unit
    @param Munit is the physical mass in units of mp per system units
    @param BoxSize half the size of a side in system units
    @param ncell is the desired number of particles per cell
    @param cnum non-zero adjusts cross section for cnum collisions per cel
    @param diamfac is a multiplicative factor for the Bohr radius (converted to system units
    @param collfrac is a multiplicative fraction for collsion number
    @param frontier enables printing of active cell list (for debugging)
    @param mfpstat enables printing of mean free path and time step statistics (for debugging)
    @param use_temp is the particle attribute position for temperature
    @param use_dens is the particle attribute position for density
    @param use_multi attempts to set the time step (if multistep is enabled)
*/
class UserTreeDSMC : public ExternalForce
{
private:
  
  int ncell;
  int cnum;
  int nsteps;
  int use_temp, use_dens;
  double diamfac;
  double diam;
  double collfrac;
  double boxsize;

  bool frontier, mfpstat, use_multi;

				// For MFP and time step statistics
  vector<double> quant, mfp_, ts_;

  unsigned stepnum;
  double curtime;

  CollideLTE *collide;
  double ElostTot;

  string comp_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg) 
  { return (NULL); }
  void initialize();

				// The tree
  double volume;

				// Timers
  Timer partnTime, treeTime;
  Timer collideTime, densityTime, gatherTime;

  void userinfo();

				// Threading
  thrd_pass_dsmc *td;
  pthread_t *t;
  void dsmc_thread_fork(pHOT *tree);


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

  //! For threading
  void *timestep_thread(void *arg);

};

#endif
