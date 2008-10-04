 #ifndef _ComponentContainer_H
 #define _ComponentContainer_H

 #include <list>
 #include <Component.H>

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
 class Interaction {
  public:
   //! Fiducial component with whom force interactions are desired
   Component* c;
   //! List of components whose particles will feel the force from *c
   list<Component*> l;
 };

//! A container for all phase-space components
/**
    There are two parameters affecting ComponentContainter, which
    may be specified in the global stanza:

    @param restart set to 1 instructs all the components and force
    methods to restart the simulation from phase space checkpoints and
    log files

    @param ctrlev is used to capture the center information for any EJ
    computations in individual components at levels smaller than the
    specifield value (default: 100)
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

public:
  //! Total number of bodies
  int ntot;

  //! Number of components
  int ncomp;

  //! List of components
  list<Component*> components;

  //! List of interactions
  list<Interaction*> interaction;

  //! Processor rates
  vector<double> rates;

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

  //! Swap states or multilevel time step interpolation
  void multistep_swap(unsigned M);

  //! Reset data before a multistep step
  void multistep_reset();

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

  //! Compute duty for each processor and initiate load balancing
  void load_balance();

};

#endif
