
#include <math.h>
#include <numerical.h>
#include <Vector.h>
#include <models.h>



double Needle::potential(Three_Vector &x)
{
	double Sp, Sm, p;

	p = b*b + SQR(x[2]) + SQR(x[3]);
	Sp = sqrt(p + SQR(x[1]+a));
	Sm = sqrt(p + SQR(x[1]-a));

	if (a < TINY)
	{
		return -G*M/sqrt(p + SQR(x[1]));
	}

	return k*log((x[1]-a+Sm)/(x[1]+a+Sp));
}


Three_Vector Needle::force(Three_Vector &x)
{
	static Three_Vector f;
	double Sp, Sm, p, q, Sprod, Ssum;

	p = b*b + SQR(x[2]) + SQR(x[3]);
	Sp = sqrt(p + SQR(x[1]+a));
	Sm = sqrt(p + SQR(x[1]-a));
	Sprod = Sp*Sm;
	Ssum = Sp+Sm;
	q = kf/p/Sprod * (Ssum - 4.0*x[1]*x[1]/Ssum);


	f[1] = -4.0*kf*x[1]/Sprod/Ssum;
	f[2] = -q*x[2];
	f[3] = -q*x[3];

	return f;
}

	

double Miyamoto_Needle::potential(Three_Vector &x)
{
	double Tp, Tm, p, q;

	p = b + sqrt(c*c + x[3]*x[3]);
	q = SQR(x[2]) + p*p;
	Tp = sqrt(q + SQR(x[1]+a));
	Tm = sqrt(q + SQR(x[1]-a));

	if (a < TINY)
	{
		return -G*M/sqrt(q + SQR(x[1]));
	}

	return k*log((x[1]-a+Tm)/(x[1]+a+Tp));
}


Three_Vector Miyamoto_Needle::force(Three_Vector &x)
{
	static Three_Vector f;
	double Tp, Tm, p, q, w, Tprod, Tsum;

	p = b + sqrt(c*c + x[3]*x[3]);
	q = SQR(x[2]) + p*p;
	Tp = sqrt(q + SQR(x[1]+a));
	Tm = sqrt(q + SQR(x[1]-a));
	Tprod = Tp*Tm;
	Tsum = Tp+Tm;
	w = kf/q/Tprod * (Tsum - 4.0*x[1]*x[1]/Tsum);


	f[1] = -4.0*kf*x[1]/Tprod/Tsum;
	f[2] = -w*x[2];
	f[3] = -w*x[3]*p/(p-b);

	return f;
}

	


double Miyamoto_Needle::density(Three_Vector &x)
{
	double ap, am, B2;


	B2 = b*b+x[2]*x[2];
	ap = a + x[1];
	am = a - x[1];

	if (a>0.0) 
	{
		return M*b/4.0/Pi/a/B2 * (
			ap/sqrt(ap*ap + B2) + am/sqrt(am*am + B2)); 
	}
	else
	{
		return M*b/2.0/Pi/pow((B2+ap*ap), 1.5);
	}
}



















