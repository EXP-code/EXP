#pragma interface

class Phase;


class Perfect
{
protected:
	double a, b;
	double alpha, beta;
	double G, M, rho0;
	double W0;
public:
	void set_Perfect(double aa, double bb, double g, double m)
	{
		a = aa;
		b = bb;
		G = g;
		M = m;
		W0 = 2.0*G*M/a/Pi;
		alpha = -a*a;
		beta = -b*b;
	}

	Three_Vector force(Three_Vector &);
	double potential(Three_Vector &);

	void Integrals(Phase &, double *, double *, double *);
	void meridional_force(double, double, double *, double *);
	void to_ellipsoidal(double, double, double *, double *);
	double second_integral(double, double, double, double, double);

	double Feff_mu(double, double);
	double Feff_lambda(double, double);

	double Fmu(double);
	double dFmu(double);
	double Flambda(double);
	double dFlambda(double);
};






