// This may look like C code, but it is really -*- C++ -*-

#ifndef _BarForcing_H
#define _BarForcing_H

#include "Perturbation.H"

/**
   Compute quadrupole expansion coefficients for a rotating ellipsoid bar
*/
class BarForcing : public Perturbation 
{

private:

  double mass, amp, length, corot, inertia;
  bool user_omega;

public:

  //! "l" harmonic index for perturbation
  static int L0;
  //! "m" harmonic index for perturbation
  static int M0;
  //! Mass of bar
  static double AMPLITUDE;
  //! Length of bar
  static double LENGTH;
  //! Minimum and maximum radii for expansion
  //@{
  static double RMIN, RMAX;
  //@}
  
  //! Constructor
  BarForcing(int NMAX, double mass, double length, double corot=1.0,
	     double amp=1.0);

  //! Parameter setup
  void compute_quad_parameters(double a21=0.5, double a32=0.1);

  //! Override omega computation
  void set_omega(double Omega) { omega = Omega; user_omega=true; }

  //! Make quadrupole bar coefficients
  void compute_coefficients();

  //! Set frequency
  void compute_omega();

  //! Return moment of intetia
  double get_Iz() { return inertia; }

  //! Return potential evaluation
  double eval(double);
};

#endif /* _BarForcing_h */
