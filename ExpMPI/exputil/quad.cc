
#include <math.h>
#include <numerical.h>
#include <Vector.h>
#include <models.h>


double Quadrupole_Bar::potential(Three_Vector &x)
{
	double s2;

	s2 = a2 + x*x;

	return -0.5*GQ*(3.0*x[1]*x[1] - x*x)/s2/s2/sqrt(s2);
}


Three_Vector Quadrupole_Bar::force(Three_Vector &x)
{
	static Three_Vector f;
	double r2, s2, s5, s7;

	r2 = x*x;
	s2 = r2 + a2;
	s5 = s2*s2*sqrt(s2);
	s7 = s5*s2;

	f = 0.5*GQ*x/s7 * (3.0*r2 - 15.0*x[1]*x[1] - 2.0*a2);
	f[1] += 3.0*GQ*x[1]/s5;

	return f;
}

	

