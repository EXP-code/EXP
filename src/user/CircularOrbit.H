// This may look like C code, but it is really -*- C++ -*-

#ifndef _CIRCULARORBIT_H

#define _CIRCULARORBIT_H 1

#include <Perturbation.H>

//! Perturbation by mass on a circular orbit
class CircularOrbit : public Perturbation {

private:

  double mass, radius;

public:
  //! Strength of perturbation
  static double AMPLITUDE;

  //! Constructor
  CircularOrbit(int NMAX, int L, int M, double mass, double radius);

  //! Compute expansion
  void compute_coefficients();

  //! Determine pattern speed
  void compute_omega();

  //! Evaluate perturbation potential
  double eval(double r);

};

#endif /* _CIRCULARORBIT_H */
