#ifndef _ComponentContainer_H
#define _ComponentContainer_H

#include <list>
#include "Component.H"
#include "Timer.H"

//! Class to keep track of the interactions between components
/** Maintains a list of component that feel the force from 
    the given component.  Every pair of interactions must be
    explicitly be listed in the [interaction] stanza of the 
    input file.  Example:
    <BR>
    <code>
    halo : disk
    </code>
    will tell the code to apply the force from halo particles
    on the disk particles.  The fully self-consistent simulation
    with a <code>disk</code> and <code>halo</code> component
    would have the lines:
    <BR>
    <code>
    halo : disk
    disk : halo
    </code>
*/
class Interaction 
{
public:
  //! Fiducial component with whom force interactions are desired
  Component* c;

  //! List of components whose particles will feel the force from *c
  std::list<Component*> l;

  //! Constructor
  Interaction(Component *c) : c(c) {}
};

//! A container for all phase-space components
/**
   There are two parameters affecting ComponentContainter, which
   may be specified in the global stanza:
   
   @param restart set to 1 instructs all the components and force
   methods to restart the simulation from phase space checkpoints and
   log files

  <br>
  Restarting:

  If the <em>infile</em> parameter contains a valid PSP file name, the
  simulation will resume using the file as an PSP body file.  If the
  <em>infile</em> begins with or contains the string 'SPL', a
  split-PSP format will be assumed.
*/
class ComponentContainer
{
  friend class Component;

private:

  bool binary;
  bool gottapot;
  double time0;
  int ifile;
  double *gcom1, *gcov1;

  //! Counts steps at each level for timing output
  vector<unsigned> levcnt;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Performance timers (enabled with VERBOSE>3)
  //@{
  Timer timer_posn, timer_gcom, timer_angmom;
  Timer timer_zero, timer_accel, timer_inter;
  Timer timer_force, timer_fixp, timer_extrn;
  Timer timer_thr_acc, timer_thr_int, timer_thr_ext;
  Timer timer_clock, timer_expand, timer_wait, timer_cuda, timer_orient;
  vector< pair<string, Timer> > timer_sntr, timer_sext;
  //@}

  //! For timers in PotAcc to hook the type of acceleration
  enum State {NONE, SELF, INTERACTION, EXTERNAL} state;

  //! Seconds between timer dumps
  static long tinterval;

  //! Timing in use
  bool timing, thread_timing;

  //! Total number of bodies
  int ntot;

  //! Number of components
  int ncomp;

  //! List of components
  std::list<Component*> components;

#ifdef HAVE_LIBCUDA
  //! Cuda particle fetch map; keep track of device particles that
  //! have been copied to host to prevent cudaMemcpy duplication
  std::map<Component*, bool> fetched;
#endif

  //! List of interactions
  std::list<std::shared_ptr<Interaction>> interaction;

  //! Processor rates
  std::vector<double> rates;

  //! Constructor
  ComponentContainer();

  //! Destructor
  ~ComponentContainer();

  //! Read database and construct indvidual components
  void initialize();

  //! Compute the acceleration for all components in the container
  void compute_potential(unsigned mlevel=0);

  /** Compute the expansion for all components in the container
      Those without expansions, should return without working */
  void compute_expansion(unsigned mlevel=0);

  //! Reset data before a multistep step
  void multistep_reset();

  /** Print out the level lists for all components to stdout 
      for diagnostic purposes */
  void print_level_lists(double T);

  //! Print column explanation for level lists output
  void print_level_list_header();

  //! Print debug info
  void multistep_debug();

  //! Debug for bad phase space values
  bool bad_values();

  /**
     Recompute the center of mass and centers of velocity for all
     components.  These may be used (if flagged) to center the phase
     space on each grid either globally or per component */
  void fix_positions();

  //! Remove the center of mass accelration from the overall phase space
  void fix_acceleration();

  //! Read initial processor rates
  void read_rates();

  //! Diagnostic report of number from each component and process
  void report_numbers();

  //! Compute duty for each processor and initiate load balancing
  void load_balance();

};

#endif
