/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the exponential integral E_1(x)
 *  (e.g. Abromiwitz and Stegun 5.1.1) using A&S 5.1.53 and 5.1.56
 *
 *
 *  Call sequence:
 *  -------------
 *  double ei1(x)
 *
 *  #include <math.h>
 *  double x;
 *
 *  Parameters:
 *  ----------
 *
 *  x        as above
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 3/4/92
 *
 ***************************************************************************/
#include <math.h>

#include "cutil.h"

#undef MAIN

double ei1(double);

static double coef0[6] = {
  -0.57721566,
   0.99999193,
  -0.24991055,
   0.05519968,
  -0.00976004,
   0.00107857 };

static double coef1a[5] = {
   0.2677737343,
   8.6347608925,
  18.0590169730,
   8.5733287401,
   1.0000000000  };

static double coef1b[5] = {
   3.9584969228,
  21.0996530827,
  25.6329561486,
   9.5733223454,
   1.0000000000  };



#ifdef MAIN
#include <stdio.h>

void usage(void){};
double series(double);
double factrl(int);

main(void)
{
  int i,numx=500;
  double xmin=0.01,xmax=20.0,dx,x;

  dx = (xmax-xmin)/(double)(numx-1);
  for (i=0; i<500; i++) {
    x = xmin + dx*i;
    printf("%e %e %e %e\n",x,ei1(x),series(x),exp(-x)/x);
  }

}

double series(x)
double x;
{
  int i,n=5,sgn=1;
  double ans=-0.5772156649;

  for (i=1; i<=n; i++) {
    sgn *= -1;
    ans -= sgn*pow(x,(double)i)/(factrl(n)*n);
  }

  ans -= log(x);
  if (ans<0.0) ans = 1.0e-20;
  return ans;
}

double factrl(n)
int n;
{
  if (n<=1)
    return 1;
  else
    return n*factrl(n-1);
}

#endif /* MAIN */

double ei1(double x)
{
  int i;
  double ans1 = 0.0;
  double ans2 = 0.0;
  
  if (x<=0.0) myerror ("ei1: x must be > 0!\n");

  if (x<=1.0) {
    for (i=5; i>0; i--) ans1 = x*(coef0[i] + ans1);
    ans1 += coef0[0];
    return ans1 - log(x);
  }

  for (i=4; i>0; i--) ans1 = x*(coef1a[i] + ans1);
  ans1 += coef1a[0];
  for (i=4; i>0; i--) ans2 = x*(coef1b[i] + ans2);
  ans2 += coef1b[0];

  return exp(-x)*ans1/(ans2*x);
}
