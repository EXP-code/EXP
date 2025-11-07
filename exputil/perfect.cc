
#include <math.h>
#include <stdio.h>
#include "numerical.H"

#include "phase.H"
#include "staeckel.H"
#include "models.H"


/*
  Potential functions for perfect ellipsoid
  
  Modified 9/24/91 to have unit mass.
*/


/*
  I always call fabs() before sqrt() here, because I had trouble with
  roundoff error as tau+alpha went to zero.
*/

double Perfect::Fmu(double mu)
{
  double x;

  x = sqrt(fabs((mu + alpha)/alpha));
  if (x <= TINY) return W0;
  
  return W0 * atanh(x)/x;
}

double Perfect::Flambda(double lambda)
{
  double x;

  x = sqrt(fabs((lambda + alpha)/alpha));
  if (x <= TINY) return W0;
  
  return W0 * atan(x)/x;
}

double Perfect::dFmu(double mu)
{
  double x, dx;
  
  x = sqrt(fabs((mu + alpha)/alpha));
  if (x <= TINY) return W0/alpha/3.0;
  dx = 0.5/x/alpha;
  
  return W0 * (-atanh(x)/x + 1.0/(1.0 - x*x)) * dx / x;
}

double Perfect::dFlambda(double lambda)
{
  double x, dx;
  
  x = sqrt(fabs((lambda + alpha)/alpha));
  if (x <= TINY) return W0/3.0/alpha;
  dx = -0.5/x/alpha;
  
  return W0 * (-atan(x)/x + 1.0/(1.0 + x*x)) * dx / x;
}




