#include <iostream>
#include <iomanip>
#include <cmath>
#include <numerical.H>

#define FACTOR 1.6
#define NTRY 50

int zbrac(std::function<double(double)> func, double *x1, double *x2)
{
  int j;
  double f1,f2;
  
  if (*x1 == *x2) {
    std::cerr << "Bad initial range in ZBRAC" << std::endl;
    exit(-1);
  }
  
  f1=func(*x1);
  f2=func(*x2);
  for (j=1;j<=NTRY;j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabs(f1) < fabs(f2))
      f1=func(*x1 += FACTOR*(*x1-*x2));
    else
      f2=func(*x2 += FACTOR*(*x2-*x1));
  }
  return 0;
}

#undef FACTOR
#undef NTRY

