// This may look like C code, but it is really -*- C++ -*-

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include <CVector.h>
#include <clinalg.h>

#include <biorth.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>
#include <gaussQ.h>

#include <CircularOrbit.H>

extern double Ylm01(int, int);

double CircularOrbit::AMPLITUDE;

CircularOrbit::CircularOrbit(int NMAX, int L, int M,
			     double Mass, double Radius) :
  Perturbation(NMAX)
{
  mass = Mass;
  radius = Radius;
  l = L;
  m = M;
  AMPLITUDE = -4.0*M_PI*mass*Ylm01(l, m)/radius;
}


void CircularOrbit::compute_coefficients()
{
  bcoef.setsize(1, nmax);
  for (int n=1; n<=nmax; n++)
    bcoef[n] = biorth->potl(n, l, biorth->r_to_rb(radius));
  bcoef *= -4.0*M_PI*mass*Ylm01(l, m);
}

void CircularOrbit::compute_omega()
{
				// Get frequency
  omega = sqrt(model->get_dpot(radius)/radius);

  omega_computed = true;
}

double CircularOrbit::eval(double r)
{
  if (r>radius)
    return AMPLITUDE*pow(r/radius, -(1.0+l));
  else
    return AMPLITUDE*pow(r/radius, l);
}
