
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

   @param mn_model uses a Miyamoto-Nagi model for the distributed gravity.
   E.g. a disk model particle.  Default: false.

   @param a is the scale length of the Miyamoto-Nagai model for the distributed gravity. Default: 0.01

   @param b is the scale height of the Miyamoto-Nagai model for the distributed gravity. Default: 0.002

   @param pm_model uses a spherical model for the distributed gravity.
   E.g. a substructure halo profile.  Default: false.

   @param pmmodel_file is the file containing the substructure halo
   profile (default: SLGridSph.model)

   @param type is the softening type. Current types: Plummer or
   Spline.  Default type is Spline.
*/

/* provide an extended spherical model for point mass */

#include <utility>
#include <memory>
#include <vector>
#include <string>
#include <set>

#include "AxisymmetricBasis.H"
#include "GravKernel.H"
#include "massmodel.H"


class Direct : public PotAccel
{
  //! Pull in base-class overrides
  using PotAccel::determine_coefficients;

private:

  int soft_indx;

  int to_proc, from_proc;
  int max_bodies;
  int ninteract;
  int ndim;

  double *tmp_buffer, *bod_buffer;

  double soft;
  bool fixed_soft;

  //! Use Miyamoto-Nagai disk
  bool mn_model;
  
  //! MN disk scale length
  double a;
  
  //! MN disk scale height
  double b;			

  //! Use point-mass profile
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

  //! Smoothing kernel isntance
  std::shared_ptr<SoftKernel> kernel;

  //! Separations smaller than this are assumed to be zero (same particle)
  const double rtol = 1.0e-16;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! The constructor
  //! \param c0 is the instantiating caller (a Component)
  //! \param conf passes in any explicit parameters
  Direct(Component* c0, const YAML::Node& conf);

  //! The constructor
  virtual ~Direct();

  //! The main force call
  void get_acceleration_and_potential(Component*);

};

#endif

