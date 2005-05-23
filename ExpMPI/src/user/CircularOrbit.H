// This may look like C code, but it is really -*- C++ -*-

#ifndef _CIRCULARORBIT_H

#define _CIRCULARORBIT_H 1

#include <Perturbation.H>

class CircularOrbit : public Perturbation {

private:

  double mass, radius;

public:
  static double AMPLITUDE;
				// Constructor

  CircularOrbit(int NMAX, int L, int M, double mass, double radius);

  void compute_coefficients();
  void compute_omega();
  double eval(double r);

};

#endif /* _CIRCULARORBIT_H */
