#include <math.h>
#include <stdio.h>
#include "cutil.h"

/*
*	Test Spline integration routine
*/
/*

main()
{
  double	a,b,func(),p,exact(),*dvector(),splsum();
  double        xx,dx,*x,*y,*z,*y2,h;
  int		n,i;
  void          splsum2();
  
  while (1)
    {
      printf("a,b,n: ");
      if (scanf("%lf %lf %d",&a,&b,&n) != 3)
	{
	  printf("wrong number of arguments\n\n");
	  break;
	}

      x = dvector(1,n);
      y = dvector(1,n);
      z = dvector(1,n);
      y2= dvector(1,n);
      dx = (b-a)/(n-1);
      for (i=1,  xx=a; i<=n; i++, xx+=dx) {
	x[i] = xx;
	y[i] = func(xx);
      }

      p = splsum(x,y,n);
      printf("ans = %35.20e   exact = %35.20e\n",p,exact(a,b));

      splsum2(x,y,z,n);
      for (i=1,  xx=a; i<=n; i++, xx+=dx)
	printf(" %le %le\n",x[i],z[i]);

    }
  
}

double	func(x)
     double	x;
{
  return(x*x*x + x*x + x + 5.0);
}
#define  integral(z) (0.25*z*z*z*z + z*z*z/3.0 + 0.5*z*z + 5.0*z)
  double	exact(a,b)
double	a,b;
{
  return(integral(b) - integral(a));
}
 */
/*   END OF TEST ROUTINE */

/*
 *	Computes integral using spline fit
 *	with a supplied integrand
 *
 *	Synopsis:
 *		ans = splsum(x,y,n);
 *		double	x[],y[];
 *              int     n;
 *
 * 	Use:
 * 		x       domain
 * 		y       range
 * 		n       number of values y[1]-->y[n]
 * 
 *      Notes:
 *              If (n<4), trapezoidal rule is used
 *
 */

double	splsum(double *x, double *y, int n)
{
  double *y2,p,h;
  int l;

  if (n < 2) {
    fprintf(stderr, "splsum: error, can't do intgral with one grid point!\n");
    return 0.0;
  }
  else if (n < 3)
    return 0.5*( (x[2]-x[1])*(y[1]+y[2]) );
  else if (n < 4)
    return 0.5*( (x[2]-x[1])*(y[1]+y[2]) + (x[3]-x[2])*(y[2]+y[3]) );


  y2 = dvector(1,n);
  spline(x,y,n,-1.0e30,-1.0e30,y2);
  p = 0.0;
  for(l=1; l<n; l++) {
    h = x[l+1] - x[l];
    p = p + 0.5*(y[l] + y[l+1])*h - (y2[l] + y2[l+1])*h*h*h/24.0;
  }
  
  free_dvector(y2,1,n);

  return p;

}


void splsum2(double *x, double *y, double *z, int n)
{
  double *y2,h;
  int l;

  if (n < 2) {
    fprintf(stderr, "splsum2: error, can't do intgral with one grid point!\n");
    z[1] = 0.0;
    return;
  }
  else if (n < 4) {
    z[1] = 0.0;
    for (l=2; l<=n; l++)
      z[l] = z[l-1] + 0.5*(y[l-1]+y[l])*(x[l]-x[l-1]);
    return;
  }

  y2 = dvector(1,n);
  spline(x,y,n,-1.0e30,-1.0e30,y2);
  z[1] = 0.0;
  for(l=2; l<=n; l++) {
    h = x[l] - x[l-1];
    z[l] = z[l-1] + 0.5*(y[l-1] + y[l])*h - (y2[l-1] + y2[l])*h*h*h/24.0;
  }
  
  free_dvector(y2,1,n);

  return;

}

