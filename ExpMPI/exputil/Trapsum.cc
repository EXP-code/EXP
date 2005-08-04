#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <interp.h>
#include <values.h>
#include <unistd.h>
#include <stdlib.h>


double	Trapsum(const Vector& x, const Vector& y)
{
  int lo = x.getlow();
  int hi = x.gethigh();

  double h, p=0.0;
  for(int l=lo; l<hi; l++) {
    h = x[l+1] - x[l];
    p = p + 0.5*(y[l] + y[l+1])*h;
  }

  return p;

}


void Trapsum(const Vector& x, const Vector& y, Vector& z)
{
  int lo = x.getlow();
  int hi = x.gethigh();

  double h;

  z[lo] = 0.0;
  for(int l=lo+1; l<=hi; l++) {
    h = x[l] - x[l-1];
    z[l] = z[l-1] + 0.5*(y[l-1] + y[l])*h;
  }
}

