/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine provides knots and weights for Gaussian integration
 *
 *
 *  Call sequence:
 *  -------------
 *  void get_gknots(knots,n);      returns structure with knots and weights
 *  free_gknots(knots,n);          frees allocated storage for above
 *
 *  #include "sphere.h"            needed to define: struct GAUSS
 *  struct GAUSS *knots;
 *  int n;
 *
 *  Parameters:
 *  ----------
 *
 *  *knots        pointer to structure which will contain knots and weights
 *  n             order of Gaussian integration
 *
 *  Returns:
 *  -------
 *
 *  None
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *
 ***************************************************************************/

#include <math.h>
#include "cutil.h"

struct GAUSS {
  int n;
  double *x,*w;
};
void get_gknots(struct GAUSS *knots, int n);
void free_gknots(struct GAUSS *knots, int n);

/*  debug */
/*
#include <stdio.h>
main()
{
  struct GAUSS rw;
  int i,n;
  double a,b,bpa,bma,*y,accum,exact(),func();

  printf("Order: ");
  scanf("%d",&n);

  get_gknots(&rw,n);

  for (i=1; i<=rw.n; i++)
    printf(" %d:  %le  %le\n",i,rw.x[i],rw.w[i]);

  printf("Will integrate x^4 - 3.5*x^3 + 2.1*x^2 + x - 5.6 . . .\n");
  printf("Enter limits a,b: ");
  scanf("%lf %lf",&a,&b);
  
  bma = 0.5*(b-a);
  bpa = 0.5*(b+a);

  y = dvector(1,rw.n);
  for(i=1; i<=rw.n; i++)
    y[i] = func(bpa + bma*rw.x[i]);
  for(i=1,accum=0.0; i<=rw.n; i++)
    accum += rw.w[i]*y[i];
  accum *= bma;

  printf("Gauss=%le   Exact=%le\n",accum,exact(b)-exact(a));

  free_gknots(&rw,n);

}

double func(x)
double x;
{
  return (-5.6 + x*(1.0 + x*(2.1 + x*(-3.5 + x))));
}
double exact(x)
double x;
{
  return x*(-5.6/1.0 + x*(1.0/2.0 + x*(2.1/3.0 + x*(-3.5/4.0 + x/5.0))));
}
*/
/* end debug */

#ifdef GENERATE

main(nargs,args)
int     nargs;
char    *args[];
{
  struct GAUSS rw;
  int i,num;
  void get_gknots();

  if (nargs==2) {
    num = atoi(args[1]);
  }
  else {
    printf("gau: usage: gau number > my.file\n");
    exit(0);
  }

  get_gknots(&rw,num);
  
  printf("%d\n",num);
  for (i=1; i<=num; i++)
    printf("%25.16e %25.16e\n",rw.x[i],rw.w[i]);

}

#endif /* GENERATE */



#define EPS 3.0e-11

void get_gknots(struct GAUSS *knots, int n)
{
  int m,j,i;
  double z1,z,pp,p3,p2,p1;
  
  knots->n = n;
  knots->x = dvector(1,n);
  knots->w = dvector(1,n);

  m=(n+1)/2;
  for (i=1;i<=m;i++)  {
    z=cos(3.14159265358979323846*(i-0.25)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);

    knots->x[i] = -z;
    knots->x[n+1-i] = z;
    knots->w[i] = 2.0/((1.0-z*z)*pp*pp);
    knots->w[n+1-i] = knots->w[i];
  }
}

#undef EPS

void free_gknots(struct GAUSS *knots, int n)
{
  free_dvector(knots->x,1,knots->n);
  free_dvector(knots->w,1,knots->n);
  knots->n = 0;
}
