
#include <math.h>
#include <numerical.h>
#include <Vector.h>
#include <models.h>


/* modified Hubble model */


double Modified_Hubble::potential(Three_Vector &x)
{
	double r, s;
	
	r = sqrt(x*x);
	s = r/rs;
	if (r < TINY)
	{
		return -Gravity*Ms/rs;
	}

	return -Gravity*Ms/r * log(s + sqrt(s*s + 1.0));
}

Three_Vector Modified_Hubble::force(Three_Vector &x)
{
	static Three_Vector f;
	double r, s, fr;
	
	r = sqrt(x*x);
	s = r/rs;
	if (r < TINY)
	{
		f.zero();
		return f;
	}

	fr = Gravity*Ms/r/r * (r / sqrt(r*r + rs*rs) - 
		log(s + sqrt(1.0 + s*s)));
	f = x*fr/r;
	return f;
}

	
	
	













