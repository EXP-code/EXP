#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <kevin_complex.h>

#ifdef CDEBUG
int Complex::nlive = 0;
#endif

Complex conjg(const Complex &z)
{
	static Complex tmp;
	tmp.r = z.r;
	tmp.i = -z.i;
	return tmp;
}

double arg(const Complex &z)
{
	double x = z.r, y = z.i;
	static double pi = 3.14159265358979;

	if (x==0.0 && y != 0.0)
	{
		if (y>=0.0) return pi/2.0;
		else return -pi/2.0;
	}
	if (x==0.0 && y==0.0) return 0.0;
	return atan2(y, x);
}

double fabs(const Complex &z)
{
	return sqrt(z.r*z.r + z.i*z.i);
}
			
		

Complex operator+(const Complex &z1, const Complex &z2)
{
	static Complex tmp; 
	tmp.r=z1.r+z2.r; 
	tmp.i=z1.i+z2.i; 
	return tmp;
}

Complex operator-(const Complex &z1, const Complex &z2)
{
	static Complex tmp; 
	tmp.r=z1.r-z2.r; 
	tmp.i=z1.i-z2.i; 
	return tmp;
}

Complex operator*(const Complex &z1, const Complex &z2)
{
	static Complex tmp; 
	tmp.r=z1.r*z2.r-z1.i*z2.i;
	tmp.i = z1.r * z2.i + z1.i * z2.r;
	return tmp;
}

Complex operator/(const Complex &z1, const Complex &z2)
{
	static Complex tmp;
	double a;
	a = z2.r*z2.r + z2.i*z2.i;
	tmp.r = (z1.r*z2.r + z1.i*z2.i)/a;
	tmp.i = (z1.i*z2.r - z1.r*z2.i)/a;
	return tmp;
}
		
Complex operator+(const Complex &z, const double x)
{
	static Complex tmp;
	tmp.r = z.r + x;
	tmp.i = z.i;
	return tmp;
}
		

Complex operator+(const double x, const Complex &z)
{
	static Complex tmp;
	tmp.r = z.r + x;
	tmp.i = z.i;
	return tmp;
}
Complex operator-(const Complex &z, const double x)
{
	static Complex tmp;
	tmp.r = z.r - x;
	tmp.i = z.i;
	return tmp;
}
		
		
Complex operator-(const double x, const Complex &z)
{
	static Complex tmp;
	tmp.r = x - z.r;
	tmp.i = -z.i;
	return tmp;
}

Complex operator*(const Complex &z, const double x)
{
	static Complex tmp;
	tmp.r = z.r*x;
	tmp.i = z.i*x;
	return tmp;
}
Complex operator*(const double x, const Complex &z)
{
	static Complex tmp;
	tmp.r = z.r*x;
	tmp.i = z.i*x;
	return tmp;
}
			
Complex operator/(const double x, const Complex &z)
{
	static Complex tmp;
	double a;
	a = z.r*z.r + z.i*z.i;
	tmp.r = x*z.r/a;
	tmp.i = -x*z.i/a;
	return tmp;
}
			
Complex operator/(const Complex &z, const double x)
{
	static Complex tmp;
	tmp.r = z.r/x;
	tmp.i = z.i/x;
	return tmp;
}

Complex exp(const Complex &z)
{
	static Complex tmp;
	double e = exp(z.r);
	tmp.r = e * cos(z.i);
	tmp.i = e * sin(z.i);
	return tmp;
}
Complex log(const Complex &z)
{
	static Complex tmp;
	double a=z.r*z.r + z.i*z.i;

	if (a==0.0) 
	{
		puts("domain error in Complex log");
		exit(0);
	}
	tmp.r = 0.5*log(a);
	tmp.i = arg(z);
	return tmp;
}
Complex pow(const Complex &base, const Complex &e)
{
	static Complex tmp;
	tmp = exp(e*log(base));
	return tmp;
}
Complex pow(const Complex &base, const double e)
{
	static Complex tmp;
	tmp = exp(e*log(base));
	return tmp;
}
Complex pow(const double base, const Complex &e)
{
	static Complex tmp;
	tmp = exp(e*log(base));
	return tmp;
}
Complex sqrt(const Complex &z)
{
	static Complex tmp;

	if (fabs(z)==0.0) {tmp.r=0.0; tmp.i=0.0;}
	else tmp = exp(0.5*log(z));

	return tmp;
}


Complex cosh(const Complex &z)
{
	static Complex tmp;
	tmp.r = cosh(z.r)*cos(z.i);
	tmp.i = sinh(z.r)*sin(z.i);
	return tmp;
}
Complex sinh(const Complex &z)
{
	static Complex tmp;
	tmp.r = sinh(z.r)*cos(z.i);
	tmp.i = cosh(z.r)*sin(z.i);
	return tmp;
}
Complex tanh(const Complex &z)
{
	static Complex c;
	c = cosh(z);
	if (fabs(c)==0.0) 
	{
		puts("domain error in tanh");
		exit(0);
	}
	return sinh(z)/c;
}

Complex asinh(const Complex &z)
{
	static Complex tmp;
	tmp = log(z + sqrt(1.0 + z*z));
	return tmp;
}

Complex acosh(const Complex &z)
{
	static Complex tmp;
	tmp = log(z + sqrt(z*z - 1.0));
	return tmp;
}

Complex atanh(const Complex &z)
{
	if (fabs(z-1.0) == 0.0 || fabs(z+1.0) == 0.0)
	{
		puts("domain error in Complex atanh()");
		exit(0);
	}
	return 0.5*log((1.0 + z)/(1.0 - z));
}


Complex cos(const Complex &z)
{
	static Complex tmp;
	tmp.r = cos(z.r)*cosh(z.i);
	tmp.i = -sin(z.r)*sinh(z.i);
	return tmp;
}


Complex sin(const Complex &z)
{
	static Complex tmp;
	tmp.r = sin(z.r)*cosh(z.i);
	tmp.i = cos(z.r)*sinh(z.i);
	return tmp;
}

Complex tan(const Complex &z)
{
	static Complex c;
	
	c = cos(z);
	if (fabs(c)==0.0) 
	{
		puts("domain error in tan");
		exit(0);
	}
	return sin(z)/c;
}


Complex acos(const Complex &z)
{
	static Complex i(0., 1.);

	return acosh(z)/i;
}

Complex asin(const Complex &z)
{
	static Complex i(0., 1.);

	return asinh(i*z)/i;
}

Complex atan(const Complex &z)
{
	static Complex i(0., 1.);

	return atanh(i*z)/i;
}

/* */

ostream& operator << (ostream& s, Complex& x)
{
  return s << "(" << x.real() << ", " << x.imag() << ")" ;
}

istream& operator >> (istream& s, Complex& x)
{

#ifdef __GNUC__
  if (!s.ipfx(0))
#else
  if (!std::istream::sentry(s))
#endif
  {
    s.clear(ios::failbit|s.rdstate());
    return s;
  }

  double r, i;
  char ch;
  s >> ws;
  s.get(ch);
  if (ch == '(')
  {
    s >> r;
    s >> ws;
    s.get(ch);
    if (ch == ',')
    {
      s >> i;
      s >> ws;
      s .get(ch);
    }
    else
      i = 0;
    if (ch != ')')
      s.clear(ios::failbit);
  }
  else
  {
    s.putback(ch);
    s >> r;
    i = 0;
  }
  x = Complex(r, i);
  return s;
}

/* */
