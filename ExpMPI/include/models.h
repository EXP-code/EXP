#pragma interface

class Miyamoto;
class Isothermal;
class Hernquist;
class Three_Vector;
class Needle;
class Miyamoto_Needle;
class Modified_Hubble;
class Quadrupole_Bar;




/*
	class definitions for galactic components 
*/

class Miyamoto
{
protected:
	double a;
	double b;
	double GM;
	double G;
	double M;

public:
	void set_Miyamoto(double ainit, double binit, double Ginit, double Minit)
	{
		a = ainit;
		b = binit;
		G = Ginit;
		M = Minit;
		GM = G*M;
	}
	void set_a(double aa) {a=aa;}
	void set_b(double bb) {b=bb;}
	void set_M(double m) {M=m; GM=G*M;}

	double get_a(void) {return a;}
	double get_b(void) {return b;}
	double get_M(void) {return M;}


	Three_Vector force(Three_Vector &);
	double potential(Three_Vector &);
};


class Isothermal
{
protected:
	double r0;
	double vc;
	double vcsq;
public:
	void set_Isothermal(double r0init, double vcinit)
	{
		r0 = r0init;
		vc = vcinit;
		vcsq = vc*vc;
	}

	Three_Vector force(Three_Vector &);
	double potential(Three_Vector &);
};


class Hernquist
{
protected:
	double r0;
	double M;
	double G;
	double GM;
public:
	void set_Hernquist(double r0init, double Ginit, double Minit)
	{
		r0 = r0init;
		G = Ginit;
		M = Minit;
		GM = G*M;
	}

	Three_Vector force(Three_Vector &);
	double potential(Three_Vector &);
};


class Logarithmic
{
protected:
	double vc;
	double vcsq;
	double r0, qy, qz;
public:
	void set_Logarithmic(double r0init, double vcinit, double qyinit, double qzinit)
	{
		r0 = r0init;
		qy = qyinit;
		qz = qzinit;
		vc = vcinit;
		vcsq = vc*vc;
	}
	void set_vc(double v) {vc = v; vcsq=vc*vc;}
	void set_r0(double r) {r0 = r;}

	double get_vc(void)   {return vc;}
	double get_r0(void)   {return r0;}

	Three_Vector force(Three_Vector &);
	double potential(Three_Vector &);
};


class Plummer
{
protected:
	double r0, r2, r3;
	double G, M, GM;
public:
	void set_Plummer(double r, double g, double m)
	{
		G = g;
		M = m;
		GM = G*M;
		r0 = r;
		r2 = r*r;
		r3 = r*r*r;
	}

	void set_M(double m) {M=m; GM=G*M;}
	void set_r0(double r) {r0=r; r2=r*r; r3=r*r*r;}

	double get_M(void)  {return M;}
	double get_r0(void) {return r0;}

	
	Three_Vector force(Three_Vector &);
	double potential(Three_Vector &);
};
		





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
		if (a > TINY) k = G*M/2.0/a;
		kf = G*M/2.0;
	}

	void set_a(double aa) {a=aa; if (a > TINY) k = G*M/2.0/a;}
	void set_b(double bb) {b=bb;}
	void set_M(double m)  {M=m; if (a>TINY) k=G*M/2.0/a; kf=G*M/2.0;}
	
	double get_a(void) {return a;}
	double get_b(void) {return b;}
	double get_M(void) {return M;}
	
	double potential(Three_Vector &);
	Three_Vector force(Three_Vector &);
};


class Point_Quadrupole
{
protected:
	double Q;
public:
	void set_Point_Quadrupole(double q)
	{
		Q = q;
	}
	double potential(Three_Vector &);
	Three_Vector force(Three_Vector &);
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

	void set_a(double aa) {a=aa; if (a>TINY) k = G*M/2.0/a;}
	void set_b(double bb) {b=bb;}
	void set_c(double cc) {c=cc;}
	void set_M(double m)  {M=m; if (a>TINY) k=G*M/2.0/a; kf=G*M/2.0;}

	double get_a(void) {return a;}
	double get_b(void) {return b;}
	double get_c(void) {return c;}
	double get_M(void) {return M;}

	double density(Three_Vector &);
	double potential(Three_Vector &);
	Three_Vector force(Three_Vector &);
};


class Modified_Hubble
{
protected:
	double rs, Ms, Gravity;
public:
	void set_Modified_Hubble(double r, double m, double G)
	{
		rs = r;
		Ms = m;
		Gravity = G;
	}

	void set_rs(double r)		{rs = r;}
	void set_Ms(double m)		{Ms = m;}
	void set_Gravity(double g)	{Gravity = g;}
	double get_rs(void) {return rs;}
	double get_Ms(void) {return Ms;}
	double get_Gravity(void) {return Gravity;}

	double potential(Three_Vector &);
	Three_Vector force(Three_Vector &);
};





class Quadrupole_Bar
{
protected:
	double GQ, a2;
public:
	void set_Quadrupole_Bar(double G, double Q, double length)
	{
		a2 = length*length;
		GQ = G*Q;
	}

	double potential(Three_Vector &);
	Three_Vector force(Three_Vector &);
};









	

	

