#ifndef _TwoCenter_H
#define _TwoCenter_H

#include <SLGridMP2.h>

#include <CBrock.H>
#include <Sphere.H>
#include <Bessel.H>
#include <CBrockDisk.H>
#include <Hernquist.H>
#include <Cylinder.H>

class MixtureBasis;

/** Computes the potential, acceleration and density using two basis
    sets with different expansion centers

    This class provides a virtual member "mixture(double *pos)" that
    must be implemented by a derived class to provide the mixture
    function for associating particles with a force center for
    gravitational field determination.  EJcom is an example of an
    implementation using an erf-based ramp whose arguement is the
    distance from the particle to the inner center and the width is
    the distance between the centers.

    All parameters are parsed by the underlying basis instances
    for the two centers.  There is only one diagnostic parameter:

    \param nhisto enables diagnostic output of mixture distribution if set to a value larger than 1
    \param basis is the name of the force method to use for the two-center expansion.  Available types are: "bessel" (Bessel), "c_brock" (CBrock), "c_brock_disk" (CBrockDisk), "hernq" (Hernquist), "sphereSL" (Sphere), "cylinder" (Cylinder).
*/
class TwoCenter : public PotAccel
{

private:

  Basis * exp_ej;
  Basis * exp_com;

  string basis;

  void initialize(void);

  MixtureBasis *mix_ej, *mix_com;

				// Will not multithread directly
  void * determine_coefficients_thread(void * arg) {};
  void * determine_acceleration_and_potential_thread(void * arg) {};

  void reset_histo();
  void write_histo();

protected:

  void accum_histo(double value);

  double dz;
  string ohisto;
  unsigned nhisto;
  vector<double> histo;

public:

  typedef double (TwoCenter::*mixFunc)(double *p);

				// Global parameters
  //! INNER center for current component
  vector<double> inner;

  //! OUTER center for current component
  vector<double> outer;

  //! Constructor
  /*! \param line passes any parameters (rsphSL, rmin and numr) */
  TwoCenter(string& line);

  //! Destructor
  ~TwoCenter();


  //! Call Basis routines to make coefficients
  //@{
  void determine_coefficients(Component *c);
  void determine_coefficients(void);
  //@}
    
  //! Main force evaluation routine: cycle through each grid
  void get_acceleration_and_potential(Component*);

  //! Mixture function for the two expansions
  //@{
  /// To be provided: must return a value in [0,1]
  virtual double mixture(double* p) = 0;
  /// The complement of the mixture function, also in [0,1]
  double Cmixture(double *p) { return 1.0 - mixture(p); }
  //@}
  
  //! For access to parent component
  virtual void RegisterComponent(Component *c) { 
    component = c; 
    exp_ej  -> RegisterComponent(c);
    exp_com -> RegisterComponent(c);
  }

  //! Multistepping routines: call Basis members
  //@}
  virtual void set_multistep_level(unsigned n) 
  { 
    mlevel = n;
    exp_ej  -> set_multistep_level(n);
    exp_com -> set_multistep_level(n);
  }

  virtual void multistep_update_begin()
  {
    exp_ej  -> multistep_update_begin();
    exp_com -> multistep_update_begin();
  }

  virtual void multistep_update(int cur, int next, Component* c, int i, int id) 
  {
    exp_ej  -> multistep_update(cur, next, c, i, id);
    exp_com -> multistep_update(cur, next, c, i, id);
  }

  virtual void multistep_update_finish()
  {
    exp_ej  -> multistep_update_finish();
    exp_com -> multistep_update_finish();
  }

  virtual void multistep_reset()
  {
    exp_ej  -> multistep_reset();
    exp_com -> multistep_reset();
  }

  virtual void multistep_debug()
  {
    exp_ej  -> multistep_debug();
    exp_com -> multistep_debug();
  }
  //@}

};

#endif
