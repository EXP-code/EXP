#include <unistd.h>
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <numerical.h>


/*
	These functions are extensions to the standard simple math functions.
Some are defined in one or the other of Borland C or unix C; I have
written code for these for compatibility. I have also added a few functions
of my own.

Extensions to Borland, defined in UNIX:
	double asinh(double x)	- inverse hyperbolic sine.
	double acosh(double x)	- inverse hyperbolic cosine (x>=1).
	double atanh(double x)	- inverse hyperbolic tangent (abs(x) <= 1).

Extensions to UNIX, defined in Borland:
	double poly(x, n, a) 	- a[0] + a[1] x + ... + a[n] x^n.
	double pow10(double x)	- 10^x


Home-grown extensions:
	int nint(double x)		- nearest integer to x, returned as an int.
	double round(double x)	- nearest integer to x, returned as a double.
	double arctanxy(double x, double y)	- atan(y/x), with safe singularities.

	void trigs(double x, int n, double *c, double *s) 
		- computes cos(0 x) ... cos(n x), and
		sin(0 x) ... sin(n x) with a recurrence relation.
		Only two trig function calls are needed.



	Kevin Long
	7/21/91
	ANSI Version 9/30/91
*/
	


/* 
Some rounding functions not defined in standard C. 
	nint returns the integer nearest to a double.
	round returns that integer as a double.
*/


double round(double x)
{
	double d1, d2;

	d1 = floor(x)-x;
	d2 = ceil(x) - x;
	if (fabs(d1)<fabs(d2)) return floor(x);
	else return ceil(x);
}

#ifdef __TURBOC__
int nint(double x)
{
	double d1, d2;

	d1 = floor(x)-x;
	d2 = ceil(x) - x;
	if (fabs(d1)<fabs(d2)) return (int) floor(x);
	else return (int) ceil(x);
}
#endif

/* 
	A fix for the trouble with atan2 at +/- Pi/2.
	Also, the name helps me remember the order of the arguments.
*/

double arctanxy(double x, double y)
{
    if (x==0.0 && y!=0.0) 
    {
	if (y>=0.0) return Pi/2.0;
	else return -Pi/2.0;
	}
	if (x==0.0 && y==0.0) return 0.0;
	return atan2(y, x);
}

/* 
	Some Borland math library extensions not defined in UNIX.
*/
#ifndef __TURBOC__

double poly(double x, int n, double *a)
{
	double p;
	int j;

	p = a[n];
	for (j=n-1; j>=0; j--) p = p*x+a[j];
	return p;
}

double pow10(double x)
{
    return pow(10.0, x);
}

#endif


/* 
	The inverse hyperbolic function are in the UNIX math library, 
but not in Borland's. 
*/

#ifdef __TURBOC__

double asinh(double x)
{
    return log(x + sqrt(1.0+x*x));
}

double acosh(double x)
{
	if (x<1.0)
	{
	  cerr.form("Domain error in acosh(): arg1=%12.5gl\n", x);
	  exit(1);
	}
	return log(x + sqrt(-1.0 + x*x));
}

double atanh(double x)
{
	if (fabs(x)>=1.0)
	{
	  cerr.form("Domain error in atanh(): arg1=%12.5lg\n", x);
	  exit(1);
	}
	return 0.5*log((1.0+x)/(1.0-x));
}

#endif









/*
	Fast evaluation of multiple-angle trig functions,
using recurrence relations. 





*/


void trigs(double x, int n, double *ct, double *st)
{
	int i;

	
	st[0] = 0.0;
	ct[0] = 1.0;

	if (n==0) return;

	
	st[1] = sin(x);
	ct[1] = cos(x);

	if (n==1) return;

	

	/* go through the recursion */
	
	for (i=2; i<=n; i++)
	{
		ct[i] = 2.0*ct[1]*ct[i-1] - ct[i-2];
		st[i] = 2.0*ct[1]*st[i-1] - st[i-2];
	}
}
	

