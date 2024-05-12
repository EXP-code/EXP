#ifndef _TwoCenter_H
#define _TwoCenter_H

#include <SLGridMP2.H>

#include <Sphere.H>
#include <Bessel.H>
#include <Cylinder.H>

class MixtureBasis;

/** Computes the potential, acceleration and density using two basis
    sets with different expansion centers

    This class provides a virtual member "mixture(double *pos)" that
    must be implemented by a derived class to provide the mixture
    function for associating particles with a force center for
    gravitational field determination.  EJcom is an example of an
    implementation using an erf-based ramp whose argument is the
    distance from the particle to the inner center and the width is
    the distance between the centers.

    All parameters are parsed by the underlying basis instances
    for the two centers.  There is only one diagnostic parameter:

    \param nhisto enables diagnostic output of mixture distribution if
    set to a value larger than 1

    \param basis is the name of the force method to use for the
    two-center expansion.  Available types are: "bessel" (Bessel),
    "sphereSL" (Sphere), "cylinder" (Cylinder).
*/
class TwoCenter : public PotAccel
{

private:

  //@{
  // The two basis instances
  Basis * exp_in;
  Basis * exp_out;
  //@}

  //! The name of the basis (for reflection)
  string basis;

  //! Force initialization
  void initialize(void);

  //@{
  //! Mixture containers for each basis
  MixtureBasis *mix_in, *mix_out;
  //@}

  //@{
  //! The required force members, but this class does not multithread directly
  void * determine_coefficients_thread(void * arg) { return 0; }
  void * determine_acceleration_and_potential_thread(void * arg) { return 0; }
  //@}

  //! Zero out the histogram bins
  void reset_histo();

  //! Write the current histogram to the log file
  void write_histo();

protected:

  //! Make a histogram for debugging the mixture factor distribution
  void accum_histo(double value);
  //! The bin spacing
  double dz;
  //! Name for the histogram file
  string ohisto;
  //! Number of bins
  unsigned nhisto;
  //! The histogram bins
  vector<double> histo;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;


public:

  //! The mixture function
  typedef double (TwoCenter::*mixFunc)(double *p);

				// Global parameters
  //! INNER center for current component
  vector<double> inner;

  //! OUTER center for current component
  vector<double> outer;

  /** Constructor
      \param c0 is the instantiating caller (a Component)
      \param line passes any parameters (rsphSL, rmin and numr) 
  */
  TwoCenter(Component* c0, const YAML::Node& conf);

  //! Destructor
  ~TwoCenter();


  //@{
  //! Call Basis routines to make coefficients
  void determine_coefficients(Component *c);
  void determine_coefficients(void);
  //@}
    
  //! Main force evaluation routine: cycle through each grid
  void get_acceleration_and_potential(Component*);

  //@{
  //! Mixture function for the two expansions

  /// To be provided: must return a value in [0,1]
  virtual double mixture(double* p) = 0;

  /// The complement of the mixture function, also in [0,1]
  double Cmixture(double *p) { return 1.0 - mixture(p); }
  //@}
  
  //! For access to parent component
  virtual void RegisterComponent(Component *c) { 
    component = c; 
    exp_in  -> RegisterComponent(c);
    exp_out -> RegisterComponent(c);
  }

  //@{
  //! Multistepping routines: call the basis members for each component
  virtual void set_multistep_level(unsigned n) 
  { 
    mlevel = n;
    exp_in  -> set_multistep_level(n);
    exp_out -> set_multistep_level(n);
  }

  virtual void multistep_update_begin()
  {
    exp_in  -> multistep_update_begin();
    exp_out -> multistep_update_begin();
  }

  virtual void multistep_update(int cur, int next, Component* c, int i, int id) 
  {
    exp_in  -> multistep_update(cur, next, c, i, id);
    exp_out -> multistep_update(cur, next, c, i, id);
  }

  virtual void multistep_update_finish()
  {
    exp_in  -> multistep_update_finish();
    exp_out -> multistep_update_finish();
  }

  virtual void multistep_reset()
  {
    exp_in  -> multistep_reset();
    exp_out -> multistep_reset();
  }

  virtual void multistep_debug()
  {
    exp_in  -> multistep_debug();
    exp_out -> multistep_debug();
  }
  //@}

};

#endif
