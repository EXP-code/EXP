/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes Laguerre knots and weights for Gaussian integration
 *  for generalized Laguerre:
 *
 *  (alpha)
 *  L      (x) ortho. over the interval (0,infinity) with weighting function:
 *   n
 *           -x  alpha
 *   w(x) = e   x      .   Evaluation technique is based on Theorems (3.6.3),
 *  (3.6.12), (3.6.20) and (3.6.21) in Introduction to Numerical Analysis by 
 *  Stoer and Bulirsch (Springer-Verlag 1980).  The implimentatin of the QL 
 *  algorithm is from Numerical Recipes.  The normalization (p_0|p_0) is given
 *  by Lanczos' series approximation with |eps| < 10^-10.
 *
 *
 *  Call sequence:
 *  -------------
 *  roots = dvector(1,n);
 *  weights = dvector(1,n);
 *
 *  void knots(roots,weights,alpha,n);
 *
 *  #include "cutil.h"
 *  double *roots,*weights
 *
 *  double x;
 *
 *  Parameters:
 *  ----------
 *
 *  None
 *
 *  Returns:
 *  -------
 *
 *  For order n integration, pass arrays as shown and knots and weights are
 *  returned.  
 *
 *  Notes:
 *  -----
 *  Values are at least as accurate as tabulated values in A&S.
 *
 *  By:
 *  --
 *
 *  MDW 05/31/89
 *
 ***************************************************************************/


#include <stdio.h>
#include <math.h>
#include "/home/weinberg/include/cutil.h"

/* TEST ROUTINE */
/*
#define QTEST 1.3459
#define LTEST 8.3323

main()
{
  double alpha,*roots,*weights,*dvector();
  double testfunc1(),testfunc2(),testfunc3(),testfunc4(),ans1,ans2,ans3,ans4;
  int j,n;
  void knots();

  printf("Order, alpha: ");
  scanf("%d %lf",&n,&alpha);

  roots = dvector(1,n);
  weights = dvector(1,n);
  knots(roots,weights,alpha,n);

  ans1 = ans2 = ans3 = ans4 = 0.0;
  for (j=1; j<=n; j++) {
    printf("%3d>  %22.14le  %22.14le  %22.14le\n",
	   j,roots[j],weights[j],weights[j]*exp(roots[j]));
    ans1 += weights[j]*testfunc1(roots[j],alpha);
    ans2 += weights[j]*testfunc2(roots[j],alpha);
    ans3 += weights[j]*testfunc3(roots[j],alpha);
    ans4 += weights[j]*testfunc4(roots[j],alpha);
  }

  printf("Test integral 1=%22.14le\n",ans1);
  printf("          Exact=%22.14le\n",
	 (QTEST*cos(LTEST)+sin(LTEST))/(QTEST*QTEST+1.0));
  printf("Test integral 2=%22.14le\n",ans2);
  printf("          Exact=%22.14le\n",M_PI*M_PI/6.0);
  printf("Test integral 3=%22.14le\n",ans3);
  printf("          Exact=%22.14le\n",M_PI*M_PI*M_PI*M_PI*7.0/120.0);
  printf("Test integral 4=%22.14le\n",ans4);
  printf("          Exact=%22.14le\n",0.75*sqrt(M_PI));

}

double testfunc1(x,alpha)
double x,alpha;
{
  return sin(QTEST*x + LTEST)*pow(x,-alpha);
}

#undef QTEST
#undef LTEST

double testfunc2(x,alpha)
double x,alpha;
{
  return pow(x,1.0-alpha)/(1.0-exp(-x));
}

double testfunc3(x,alpha)
double x,alpha;
{
  return pow(x,3.0-alpha)/(1.0+exp(-x));
}

double testfunc4(x,alpha)
double x,alpha;
{
  return pow(x,1.5-alpha);
}

*/
/* END TEST ROUTINE */

#ifdef GENERATE

main(nargs,args)
int     nargs;
char    *args[];
{
  double *eroots,*ewghts,power,alpha;
  int i,num;
  void knots();

  if (nargs==3) {
    num = atoi(args[1]);
    alpha = atof(args[2]);
  }
  else {
    printf("leg: usage: leg number power > my.file\n");
    exit(0);
  }

  eroots = dvector(1,num);
  ewghts = dvector(1,num);
  knots(eroots,ewghts,alpha,num);
  
  printf("%d %25.16e\n",num,alpha);
  for (i=1; i<=num; i++)
    printf("%25.16e %25.16e\n",eroots[i],ewghts[i]);

}

#endif /* GENERATE */




void knots(double *roots, double *weights, double C, int n)
{
  double *d,*e,**z,norm,dgammln(double xx);
  int i,j;
  void tqli(double *d, double *e, int n, double **z);

  d = dvector(1,n);
  e = dvector(1,n);
  z = dmatrix(1,n,1,n);
  for (i=1; i<=n; i++) {
    for (j=1; j<=n; j++)
      z[i][j] = 0.0;
    z[i][i] = 1.0;
    d[i] = 2.0*i - 1.0 + C;
    e[i] = sqrt(((double)(i-1) + C)*(double)(i-1));
  }
  
  tqli(d,e,n,z);

  for (i=1; i<=n; i++) {
    norm = 0.0;
    for (j=1; j<=n; j++)
      norm += z[j][i]*z[j][i];
    weights[i] = z[1][i]*z[1][i]*exp(dgammln(1.0+C))/norm;
    roots[i] = d[i];
  }

  free_dvector(d,1,n);
  free_dvector(e,1,n);
  free_dmatrix(z,1,n,1,n);

}



#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
  
  void tqli(double *d, double *e, int n, double **z)
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  void nrerror();
  
  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if (fabs(e[m])+dd == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) nrerror("Too many iterations in TQLI");
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=sqrt((g*g)+1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  if (fabs(f) >= fabs(g)) {
	    c=g/f;
	    r=sqrt((c*c)+1.0);
	    e[i+1]=f*r;
	    c *= (s=1.0/r);
	  } else {
	    s=f/g;
	    r=sqrt((s*s)+1.0);
	    e[i+1]=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  p=s*r;
	  d[i+1]=g+p;
	  g=c*r-b;
	  /* Next loop can be omitted if eigenvectors not wanted */
	  for (k=1;k<=n;k++) {
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
	}
	d[l]=d[l]-p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}

