#include <cmath>

double poly(double x, int n, double *a)
{
  double p = a[n];
  for (int j=n-1; j>=0; j--) p = p*x+a[j];
  return p;
}

