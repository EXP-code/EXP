// This may look like C code, but it is really -*- C++ -*-

#ifndef _BarForcing_H

#define _BarForcing_H 1

#include "Perturbation.H"

class BarForcing : public Perturbation 
{

private:

  double mass, length, corot, inertia;
  bool user_omega;

public:

  static int L0;
  static int M0;
  static int Mnum;
  static double AMPLITUDE;
  static double LENGTH;
  static double RMIN, RMAX;

				// Constructor

  BarForcing(int NMAX, double mass, double length, double corot=1.0);

				// Setup

  void compute_quad_parameters(double a21=0.5, double a32=0.1);
  void set_omega(double Omega) { omega = Omega; user_omega=true; }

				// Make quadrupole bar coefficients
  void compute_coefficients();

				// Set frequency
  void compute_omega();

  double get_Iz() { return inertia; }
};

#endif /* _BarForcing_h */
