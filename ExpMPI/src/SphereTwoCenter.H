#ifndef _SphereTwoCenter_H
#define _SphereTwoCenter_H

#include <SLGridMP2.h>
#include <SphericalBasisMixtureSL.H>

/** Computes the potential, acceleration and density using the Sturm
    Liouville direct solution

    \param rmin is the minimum value in the table for the 
    radial basis functions (default is 0.001)

    \param rmax is the minimum value in the table for the 
    radial basis functions (default is 0.001)

    \param cmap should be 1 for coordinate scaling, 0 for no scaling

    \param Lmax is maximum order for the harmonic expansion (\f$Y_{lm}\f$)

    \param nmax is radial maximum order in the expansion basis

    \param rs is the halo scale length (default is 0.067*rmax)

    \param numr is the number radial coordinate knots in the table 
    (default is 1000)

    \param self_consistent = 0 flags the expansion to fix the coefficients
    after the first time step
*/
class SphereTwoCenter : public PotAccel
{

  friend class SphericalBasisMixtureSL;

private:

  SphericalBasisMixtureSL *exp_ej, *exp_com;

  void initialize(void);

				// Parameters
  double rsphSL;
  double rmax;
  double rmin;
  double rs;
  int numr;
  int cmap;
  int Lmax;
  int nmax;

				// Will call SphericalBasisMixtureSL members
				// to do this
  void determine_coefficients(void) {};
  
				// Will not multithread directly
  void * determine_coefficients_thread(void * arg) {};
  void * determine_acceleration_and_potential_thread(void * arg) {};

  //! Flag whether or not a call has taken place yet
  bool firstime_accel;

  //! If set, coefficients will only be computed once
  bool self_consistent;

protected:

  //! Points at the instance of SL grid shared between SphericalBasis instances
  SLGridSph* ortho;

public:

				// Global parameters
  //! Center for current component
  double *center;

  //! Constructor
  /*! \param line passes any parameters (rsphSL, rmin and numr) */
  SphereTwoCenter(string& line);

  //! Destructor
  ~SphereTwoCenter();

  //! Main force evaluation routine: cycle through each grid
  void get_acceleration_and_potential(vector<Particle>*);

  //! Mixture function for the two expansions
  virtual double mixture(Particle& p) = 0;

  //! For access to parent component
  virtual void RegisterComponent(Component *c) { 
    component = c; 
    exp_ej->RegisterComponent(c);
    exp_com->RegisterComponent(c);
  }

};

#endif


