
#include <stdio.h>
#include <math.h>
#include <numerical.h>
#include <Vector.h>

#include <phase.h>
#include <models.h>




Three_Vector Isothermal::force(Three_Vector &u)
{
	double r, fr, x;
	static Three_Vector tmp;

	r = sqrt(u*u);

	if (r <= TINY)
	{
		tmp.zero();
		return tmp;
	}

	x = r/r0;

	fr =  -vcsq*r0*(x - atan(x))/r/r;

	tmp = fr*u/r;

	return tmp;
}

double Isothermal::potential(Three_Vector &u)
{
	double x, r;

	r = sqrt(u*u);
	x = r/r0;

	return vcsq*(atan(x)/x + log(1.0+x*x)/2.0);
}




Three_Vector Miyamoto::force(Three_Vector &u)
{
	double rc, z, s, t, w;
	static Three_Vector tmp;

	tmp = u;

	rc = sqrt(SQR(u[1]) + SQR(u[2]));	/* cylindrical radius */
	z = u[3];

	t = sqrt(z*z + b*b);
	s = a + t;
	w = sqrt(s*s + rc*rc);
	w = w*w*w;

	tmp *= -GM/w;
	tmp[3] *= s/t;

	return tmp;
}

double Miyamoto::potential(Three_Vector &u)
{
	double rc, z, s, t, w;

	rc = sqrt(SQR(u[1]) + SQR(u[2]));	/* cylindrical radius */
	z = u[3];
	t = sqrt(z*z + b*b);
	s = a + t;
	w = sqrt(s*s + rc*rc);

	return -GM/w;
}


Three_Vector Hernquist::force(Three_Vector &u)
{
	double r, fr;
	static Three_Vector tmp;

	r = sqrt(u*u);

	fr = -GM/SQR(r+r0);

	tmp = fr*u/r;
	return tmp;
}


double Hernquist::potential(Three_Vector &u)
{
	double r;
	r = sqrt(u*u);

	return -GM/(r + r0);
}


double Logarithmic::potential(Three_Vector &u)
{
	return 0.5*vcsq*log(r0*r0 + u[1]*u[1] + SQR(u[2]/qy) 
		+ SQR(u[3]/qz));
}

Three_Vector Logarithmic::force(Three_Vector &u)
{
	double s;
	static Three_Vector tmp;
	tmp = u;

	s = r0*r0 + u[1]*u[1] + SQR(u[2]/qy) + SQR(u[3]/qz);

	tmp *= -vcsq/s;
	tmp[2] *= 1.0/qy/qy;
	tmp[3] *= 1.0/qz/qz;

	return tmp;
}



Three_Vector Plummer::force(Three_Vector &x)
{
	static Three_Vector f;
	double s;

	s = sqrt(1.0 + x*x/r2);

	f = -GM*x/s/s/s/r3;

	return f;
}

double Plummer::potential(Three_Vector &x)
{
	return -GM/r0/sqrt(1.0 + x*x/r2);
}



	
	
	
	







	
















