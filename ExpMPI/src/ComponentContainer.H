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
class ComponentContainer
{
  friend class Component;

private:

  bool binary;
  double time0;
  string infile;
  string outfile;
  int ifile;

public:
  //! Total number of bodies
  int ntot;

  //! Number of components
  int ncomp;

  //! List of components
  list<Component*> components;

  //! List of interactions
  list<Interaction*> interaction;

  //! Constructor
  ComponentContainer();

  //! Destructor
  ~ComponentContainer();

  //! Read database and construct indvidual components
  void initialize();

  //! Compute the acceleration for all components in the container
  void compute_potential();

  /**
     Recompute the center of mass and centers of velocity for all
     components.  These may be used (if flagged) to center the phase
     space on each grid either globally or per component */
  void fix_positions();

  //! Remove the center of mass accelration from the overall phase space
  void fix_acceleration();

};

#endif
