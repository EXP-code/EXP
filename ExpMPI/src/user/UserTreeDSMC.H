#ifndef _UserTreeDSMC_H
#define _UserTreeDSMC_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalForce.H>

#include <RunningTime.H>

class CollideLTE;

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
    @param tsdiag enables printing gas time step diagnostics (for debugging)
    @param mfpstat enables printing of mean free path and time step statistics (for debugging)
    @param cbadiag enables printing of length-scale statistics for the CBA algorithm (for debugging)
    @param dryrun enabled turns off collsions (primarily for debugging)
    @param nsteps is the number of steps between diagnostic output (default: -1)
    @param msteps is the maximum level for diagnostic output (default: -1)
    @param use_temp is the particle attribute position for temperature
    @param use_dens is the particle attribute position for density
    @param use_multi attempts to set the time step (if multistep is enabled)
    @param use_repair repairs the tree rather than remakes the tree
    @param use_delt include energy loss rate in setting time step (if multistep is enabled) 
    @param use_exes track energy loss per particle for exact energy conservation
    @param use_Kn track Knudsen number in attribute position n (default: -1)
    @param use_St track Strouhal number in attribute position n (default: -1)
    @param esol solves the LTE cooling rate equation for each cell
    @param coolfrac is the fraction of the cooling time for setting the time step (default: 0.1)
    @param cba use the Collisional Boltzmann Approximation from Alexander et al.
    @param madj levels and below perform a full domain decomposition and above madj the tree is repaired (default: 0)
    @param wght sets the cell time partitioning (default: 1 = on)
    @param sub_sample toggles the key partition algorithm that uses only a subsample of the key list to decompose the domaine (default: true)
*/
class UserTreeDSMC : public ExternalForce
{
private:
  
  int ncell, Ncell;
  int cnum, wght;
  int nsteps, msteps;
  unsigned madj;
  int use_temp, use_dens, use_delt, use_exes, use_Kn, use_St, tspow;
  double diamfac;
  double diam;
  double collfrac;
  double boxsize;
  double boxratio;
  double jitter;
  double epsm;
  double coolfrac;

  bool frontier, tsdiag, mfpstat, statall, use_multi, use_pullin, use_repair;
  bool cba, cbadiag, dryrun, nocool, ntc, tube, esol, sub_sample;

				// For MFP and time step statistics
  vector<double> quant, mfp_, ts_, coll_, nsel_, cool_, rate_;

  unsigned stepnum;
  unsigned llevel, clevel;
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
  Timer partnTime, tree1Time, tree2Time, tstepTime, llistTime, collideTime;
  RunningTime overhead;

  void userinfo();


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

};

#endif
