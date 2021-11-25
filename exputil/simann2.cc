// simanneal.c++   Implementation of a General Purpose Simulated Annealing Class
/* Uses Cauchy training         */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <random>
#include <cmath>

#include <simann2.h>

extern std::mt19937 random_gen;

#ifndef HUGE
#define HUGE	HUGE_VAL
#endif

SimAnneal::SimAnneal (Func1d f, const int d):
  func (f), dimension (d), ddwell (20), rrange (M_PI/2.0), t0 (0.0), K (1.0),
  rho (0.5), dt (0.1), tscale (0.1), maxit (400), c_jump (100.0), fsave(0)
{
  std::uniform_real_distribution<>::param_type params1(-rrange, rrange);
  std::uniform_real_distribution<>::param_type params2(0.0, 1.0);

  number_range.param(params1);
  number_01.param(params2);

  x    .resize(dimension);
  xnew .resize(dimension);
  xbest.resize(dimension);

  y = ybest = HUGE;

  err = 0;
}

int SimAnneal::
set_up (Func1d f, const int d, const uint32_t seed)
{

  dimension = d;

  func = f;

  x     .resize(dimension);
  xnew  .resize(dimension);
  xbest .resize(dimension);

  y = ybest = HUGE;

  err = 0;

  return err;
}

// increase the temperature until the system "melts"
double SimAnneal::
melt (const int iters)
{
  int ok = 0;
  double xc, ynew, t, cold, c = 0.0;

  int n = iters;
  if (n < 1)
    n = maxit;

  t = t0;

  for (int i=0; i<n; i++) {

    if (i > 0 && c > 0.0) {
      cold = c;
      ok = 1;
    }
    
    t += dt;

    for (int j=0; j<dimension; j++) {
      xc = rho * t * tan (number_range(random_gen));
      x[j] += xc;
    }

    equilibrate (t, ddwell);

    // "energy"
    ynew = func(x);
    c = ynew - y;
    
    if (c < 0.0 && ynew < ybest) {
      for (int j=0; j<dimension; j++)
	xbest[j] = x[j];
      ybest = ynew;
    }

    y = ynew;

    if (ok && c > (c_jump * cold))	/* phase transition */
      break;
    
  }

  return t0 = t;
}

// iterate a few times at the present temperature to get to thermal
// equilibrium
//
int SimAnneal::
equilibrate (const double t, const int n)
{
  int equil = 0;
  double xc, ynew, c, delta;

  delta = 1.0;
  int i = 0;
  for (; i<n; i++) {
    for (int j=0; j<dimension; j++) {
      xc = rho * t * tan (number_range(random_gen));
      xnew[j] = x[j] + xc;
    }

    // "energy"
    ynew = func(xnew);
    c = ynew - y;

    if (c < 0.0) {		// keep xnew if energy is reduced 
      std::swap(x, xnew);
      
      y = ynew;
      if (y < ybest) {
	for (int j=0; j<dimension; j++) xbest[j] = x[j];
	ybest = y;
      }

      delta = fabs (c);
      delta = (y != 0.0) ? delta / y : (ynew != 0.0) ?
	delta / ynew : delta;
      
      // equilibrium is defined as a 10% or smaller change in 10
      // iterations
      if (delta < 0.10)
	equil++;
      else
	equil = 0;
    }
    else {
      // keep xnew with probability, p, if ynew is increased
      /*
	p = exp( - (ynew - y) / (K * t) );

	if ( p > number_01(random_gen) )
	{
	xtmp = x;
	x = xnew;
	xnew = xtmp;
	y = ynew;
	equil = 0;
	}
	else
      */
      
      equil++;
    }

    if (equil > 9)
      break;
  }

  return i + 1;
}

// cool the system with annealing 
double SimAnneal::
anneal (const int iters)
{
  double xc, p, ynew, t, c, dt, told;

  int n = iters;
  if (n < 1)
    n = maxit;

  equilibrate (t0, 10 * ddwell);

  told = t0;
  for (int i=0; i<n; i++) {
    t = t0 / (1.0 + i * tscale);
    dt = t - told;
    told = t;
    
    equilibrate (t, ddwell);
    
    for (int j=0; j<dimension; j++) {
      xc = rho * t * tan (number_range(random_gen));
      xnew[j] = x[j] + xc;
    }

    // "energy"
    ynew = func(xnew);
    c = ynew - y;

    if (ynew <= y) {		// keep xnew if energy is reduced
      std::swap(x, xnew);
      y = ynew;
      
      if (y < ybest) {
	for (int j=0; j<dimension; j++)
	  xbest[j] = x[j];
	ybest = y;
      }
      continue;
    } else {
      // keep xnew with probability, p, if ynew is increased
      p = exp (-(ynew - y) / (K * t));
      
      if (p > number_01(random_gen)) {
	std::swap(x, xnew);
	y = ynew;
      }
    }
    
    if (fsave) log_state(i);
  }
  
  equilibrate (t, 10 * ddwell);

  t0 = t;

  return t;
}

void SimAnneal::
log_state(int i)
{
  std::ofstream out(fname.c_str(), std::ios::app);
  out << std::setw(5) << i;
  out << std::setw(15) << y;
  for (int j=0; j<dimension; j++) out << std::setw(15) << x[j];
  out << std::endl;
}
