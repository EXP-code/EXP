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
    @param Lunit is the physical scale in cm per system units (default: 300 kpc)
    @param Tunit is the physical time scale in years per system unit (default: gravitational)
    @param Munit is the physical mass in units of mp per system units (default: 1e12 Msun)
    @param boxsize half the size of a side in system units (default: 1)
    @param boxratio is the ratio of the z box height to the x or y box height (default: 1)
    @param jitter is the fraction of the box size added to the tree origin (default: 0)
    @param ncell is the desired number of particles per cell (default: 7)
    @param Ncell is the desired number of particles per sample cell (default: 64)
    @param cnum non-zero adjusts cross section for cnum collisions per cell
    @param epsm >0 enables EPSM computation for cells with  MFP/size ratio < epsm
    @param diamfac is a multiplicative factor for the Bohr radius (converted to system units
    @param collfrac is a multiplicative fraction for collsion number
    @param frontier enables printing of active cell list (for debugging)
    @param mfpstat enables printing of mean free path and time step statistics (for debugging)
    @param dryrun enabled turns off collsions (primarily for debugging)
    @param use_temp is the particle attribute position for temperature
    @param use_dens is the particle attribute position for density
    @param use_multi attempts to set the time step (if multistep is enabled)
    @param use_delt include energy loss rate in setting time step (if multistep is enabled) 
    @param use_exes track energy loss per particle for exact energy conservation
    @param coolfrac is the fraction of the cooling time for setting the time step (default: 0.1)
    @param cba use the Collisional Boltzmann Approximation from Alexander et al.
*/
class UserTreeDSMC : public ExternalForce
{
private:
  
  int ncell, Ncell;
  int cnum;
  int nsteps;
  int use_temp, use_dens, use_delt, use_exes;
  double diamfac;
  double diam;
  double collfrac;
  double boxsize;
  double boxratio;
  double jitter;
  double epsm;
  double coolfrac;

  bool frontier, mfpstat, use_multi, use_pullin, cba, dryrun, nocool;

				// For MFP and time step statistics
  vector<double> quant, mfp_, ts_, coll_, nsel_, rate_;

  unsigned stepnum;
  double curtime;

  CollideLTE *collide;
  double ElostTotCollide, ElostTotEPSM;

  string comp_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg) 
  { return (NULL); }
  void initialize();

				// The tree
  double volume;

				// Timers
  Timer partnTime, treeTime, collideTime;

  void userinfo();

				// Threading
  thrd_pass_dsmc *td;
  pthread_t *t;
  void dsmc_thread_fork(pHOT *tree);

				// Debugging
  void triggered_cell_body_dump(double time, double radios);

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
