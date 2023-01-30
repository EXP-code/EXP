#include <math.h>
#include <numerical.H>
#include <stdio.h>


static std::function<double(double, double)> q2_fdum;
static double q2_ay, q2_by, q2_tol, q2_xdum;


/* 
	Adaptive quadrature in 2D. Uses recursive calls to 
qadapt(). See Numerical Recipes in C, pg. 140. for an analogous routine
with Gaussian quadrature.

*/

double qadapt2d(double ax, double bx, double ay, double by,
		std::function<double(double, double)> f, double tol)
{
	double f1(double);

	q2_fdum = f;

	q2_ay = ay;
	q2_by = by;
	q2_tol = tol;

	return qadapt(ax, bx, f1, tol);
}


double f1(double x)
{
	double f;
	double f2(double);

	q2_xdum = x;

	f = qadapt(q2_ay, q2_by, f2, q2_tol);
	return f;
}



double f2(double y)
{
	double f;
	f = q2_fdum(q2_xdum, y);
	return f;
}
	






