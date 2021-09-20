#ifndef _Needle_H
#define _Needle_H


class Needle
{
protected:	
	double a, b;
	double G, M, k, kf;
public:
	void set_Needle(double aa, double bb, double gg, double mm)
	{
		a = aa;
		b = bb;
		G = gg;
		M = mm;
		k = G*M/2.0/a;
		kf = G*M/2.0;
	}

	double potential(Vector &);
	Vector force(Vector &);
};


class Miyamoto_Needle
{
protected:	
	double a, b, c;
	double G, M, k, kf;
public:
	void set_Miyamoto_Needle(double aa, double bb, double cc, 
		double gg, double mm)
	{
		a = aa;
		b = bb;
		c = cc;
		G = gg;
		M = mm;
		k = G*M/2.0/a;
		kf = G*M/2.0;
	}

	double potential(Vector &);
	Vector force(Vector &);
};

#endif
