// This may look like C code, but it is really -*- C++ -*-

#pragma implementation "simann2.h"

// simanneal.c++   Implementation of a General Purpose Simulated Annealing Class
/* Uses Cauchy training         */

static const char rcsid[] = "@(#)simann.c++	1.3 15:55:47 3/30/93   EFC";

#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <simann2.h>

#ifndef HUGE
#define HUGE	HUGE_VAL
#endif


#ifndef PI2
#define PI2		(PI/2.0)
#endif


SimAnneal::SimAnneal (Func1d* f, const int d):
  func (f), dimension (d), ddwell (20), rrange (PI2), t0 (0.0), K (1.0),
  rho (0.5), dt (0.1), tscale (0.1), maxit (400), c_jump (100.0), fsave(0)
{

  gen = new ACG (10, 20);
  number_range = new Uniform (-rrange, rrange, gen);
  number_01 = new Uniform (0.0, 1.0, gen);

  x = new double[dimension];
  xnew = new double[dimension];
  xbest = new double[dimension];
  y = ybest = HUGE;

  if (x == NULL || xnew == NULL || xbest == NULL)
    err = -1;
  else
    err = 0;
}

int SimAnneal::
set_up (Func1d *f, const int d, const uint32_t seed)
{

  dimension = d;

  func = f;

  gen = new ACG (seed, 20);

  x = new double[dimension];
  xnew = new double[dimension];
  xbest = new double[dimension];
  y = ybest = HUGE;

  if (x == NULL || xnew == NULL || xbest == NULL)
    err = -1;
  else
    err = 0;

  return err;

}

/* increase the temperature until the system "melts" */
double SimAnneal::
melt (const int iters)
{
  int i, j, ok = 0;
  double xc, ynew, t, cold, c = 0.0;

  int n = iters;
  if (n < 1)
    n = maxit;

  t = t0;

  for (i = 0; i < n; i++)
    {
      if (i > 0 && c > 0.0)
	{
	  cold = c;
	  ok = 1;
	}

      t += dt;

      for (j = 0; j < dimension; j++)
	{
	  xc = rho * t * tan ((*number_range) ());
	  x[j] += xc;
	}

      equilibrate (t, ddwell);

      /* "energy" */
      ynew = func->CostFunction (x);
      c = ynew - y;

      if (c < 0.0 && ynew < ybest)
	{
	  for (j = 0; j < dimension; j++)
	    xbest[j] = x[j];
	  ybest = ynew;
	}

      y = ynew;

      if (ok && c > (c_jump * cold))	/* phase transition */
	break;

    }

  return t0 = t;

}

/* iterate a few times at the present temperature */
/* to get to thermal equilibrium */
int SimAnneal::
equilibrate (const double t, const int n)
{
  int i, j, equil = 0;
  double xc, ynew, c, delta;
  double *xtmp;
  //    double p;

  delta = 1.0;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < dimension; j++)
	{
	  xc = rho * t * tan ((*number_range) ());
	  xnew[j] = x[j] + xc;
	}
      /* "energy" */
      ynew = func->CostFunction (xnew);
      c = ynew - y;

      if (c < 0.0)		/* keep xnew if energy is reduced */
	{
	  xtmp = x;
	  x = xnew;
	  xnew = xtmp;

	  y = ynew;
	  if (y < ybest)
	    {
	      for (j = 0; j < dimension; j++)
		xbest[j] = x[j];
	      ybest = y;
	    }

	  delta = fabs (c);
	  delta = (y != 0.0) ? delta / y : (ynew != 0.0) ?
	    delta / ynew : delta;

	  /* equilibrium is defined as a 10% or smaller change
	     in 10 iterations */
	  if (delta < 0.10)
	    equil++;
	  else
	    equil = 0;


	}
      else
	{
	  /* keep xnew with probability, p, if ynew is increased */
	  /*
	     p = exp( - (ynew - y) / (K * t) );

	     if ( p > (*number_01)() )
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

/* cool the system with annealing */
double SimAnneal::
anneal (const int iters)
{
  int i, j;
  double xc, p, ynew, t, c, dt, told;
  double *xtmp;


  int n = iters;
  if (n < 1)
    n = maxit;

  equilibrate (t0, 10 * ddwell);

  told = t0;
  for (i = 0; i < n; i++)
    {
      t = t0 / (1.0 + i * tscale);
      dt = t - told;
      told = t;

      equilibrate (t, ddwell);

      for (j = 0; j < dimension; j++)
	{
	  xc = rho * t * tan ((*number_range) ());
	  xnew[j] = x[j] + xc;
	}
      /* "energy" */
      ynew = func->CostFunction (xnew);
      c = ynew - y;

      if (ynew <= y)		/* keep xnew if energy is reduced */
	{
	  xtmp = x;
	  x = xnew;
	  xnew = xtmp;
	  y = ynew;

	  if (y < ybest)
	    {
	      for (j = 0; j < dimension; j++)
		xbest[j] = x[j];
	      ybest = y;
	    }
	  continue;
	}
      else
	{
	  /* keep xnew with probability, p, if ynew is increased */
	  p = exp (-(ynew - y) / (K * t));

	  if (p > (*number_01) ())
	    {
	      xtmp = x;
	      x = xnew;
	      xnew = xtmp;
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
initial (double *xi)
{
  for (int k = 0; k < dimension; k++)
    x[k] = xi[k];
}

void SimAnneal::
current (double *xc)
{
  for (int k = 0; k < dimension; k++)
    xc[k] = x[k];
}

void SimAnneal::
optimum (double *xb)
{
  for (int k = 0; k < dimension; k++)
    xb[k] = xbest[k];
}

void SimAnneal::
log_state(int i)
{
  ofstream out(fname.c_str(), ios::app);
  out << setw(5) << i;
  out << setw(15) << y;
  for (int j = 0; j < dimension; j++)
    out << setw(15) << x[j];
  out << endl;
}
