#include <unistd.h>
#include <stdio.h>
#include <math.h>

void satellite_orbit(double T, double* X, double* Y, double* Z);

main()
{
  double TMIN, TMAX, T, X, Y, Z, dT;
  int i, num;
  FILE *out;

  if ( !(out = fopen("testsatorb.data", "w")) ) {
    fprintf(stderr, "Couldn't open output file\n");
    exit(-1);
  }

  printf("TMIN, TMAX, num: ");
  scanf("%lf %lf %d", &TMIN, &TMAX, &num);

  dT = (TMAX - TMIN)/num;

  for (i=0; i<=num; i++) {
    T = TMIN + dT*i;
    satellite_orbit(T, &X, &Y, &Z);

    fprintf(out, "%13.6e  %13.6e  %13.6e\n", X, Y, Z);
  }
  
}
