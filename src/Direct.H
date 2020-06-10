
#ifndef _Direct_H
#define _Direct_H

/**
   Computes the potential, acceleration and density using the direct summation

   Positions are passed from node to node on a ring.  Accelerations
   computed on requested particles array due to particles in component
   associated with the Direct force.  

   Settable parameters:
   @param soft_indx is the index in the particle double array for the 
   softening value.

   @param soft is a fixed softening value for all particles.  This can be 
   used even if a value is specified in the particle attribute array.

   @param pm_model uses a spherical model for the distributed gravity.
   E.g. a substructure halo profile.  Default: false.

   @param pmmodel_file is the file containing the substructure halo
   profile (default: SLGridSph.model)

   @param type is the softening type. Current types: Plummer or
   Spline.  Default type is Spline.
*/

/* provide an extended spherical model for point mass */

#include <functional>
#include <utility>

#include <AxisymmetricBasis.H>
#include <massmodel.h>

class Direct : public PotAccel
{
private:

  int soft_indx;

  int to_proc, from_proc;
  int max_bodies;
  int ninteract;
  int ndim;

  double *tmp_buffer, *bod_buffer;

  double soft;
  bool fixed_soft;

  bool pm_model;
  double diverge_rfac;
  int diverge;
  string pmmodel_file;
  SphericalModelTable *pmmodel;

  void initialize();

  void determine_coefficients(void);
  void determine_acceleration_and_potential(void);

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

  std::function<std::pair<double, double>(double, double)> sfct;

public:

  //! The constructor
  /** \param conf passes in any explicit parameters */
  Direct(const YAML::Node& conf);

  //! The constructor
  virtual ~Direct();

  //! The main force call
  void get_acceleration_and_potential(Component*);

};

#endif

