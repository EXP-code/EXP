
#ifndef _k_complex_h

#include <iostream.h>

#define _k_complex_h
#define Complex KComplex

class CVector;
class CMatrix;
class Complex;

// #define CDEBUG 1


class Complex
{
protected:
	double r;
	double i;
#ifdef CDEBUG
	static int nlive;
#endif
public:
	friend CVector;
	friend CMatrix;

	/* constructors */
	
#ifdef CDEBUG
	Complex(void) {r=0.0; i=0.0; nlive++;}
	Complex(const double x) {r=x; i=0.0; nlive++;}
	Complex(const double x, const double y) {r=x; i=y; nlive++;}
	Complex(const Complex &z) {r=z.r; i=z.i; nlive++;}
	~Complex(void) {nlive--;}
#else
	Complex(void) {r=0.0; i=0.0;}
	Complex(const double x) {r=x; i=0.0;}
	Complex(const double x, const double y) {r=x; i=y;}
	Complex(const Complex &z) {r=z.r; i=z.i;}
#endif
	

	/* access to real and imag parts */
	
	double & real(void) {return r;}
	double & imag(void) {return i;}	
	friend inline double Re(Complex z);
	friend inline double Im(Complex z);



	/* simple output */
	
	void print(void)	{cout << r << " + " << i << " I\n";
#ifdef CDEBUG
	static void count(void)	{
	  cerr << "number of live complexes: " << nlive << '\n';
	}
#endif	
	

	/* assignment */
	
	Complex &operator=(const Complex &z) {r=z.r; i=z.i; return *this;}
	Complex &operator=(const double &x) {r=x; i=0.0; return *this;}


	/* unary minus */

	Complex operator-(void)	
		{
			Complex tmp;
			tmp.r = -r;
			tmp.i = -i;
			return tmp;
		}



	/* reflexive arithmetic operators */
	
	Complex &operator+=(const Complex &z) {r+=z.r; i+=z.i; return *this;}
	Complex &operator+=(double x) {r+=x; return *this;}
	Complex &operator-=(const Complex &z) {r-=z.r; i-=z.i; return *this;}
	Complex &operator-=(double x) {r-=x; return *this;}
	Complex &operator*=(const Complex &z) 
		{
			double tmp=r; 
			r=r*z.r - i*z.i; 
			i=i*z.r + tmp*z.i; 
			return *this;
		}
	Complex &operator*=(double x) {r *= x; i *= x; return *this;}
	Complex &operator/=(double x) {r /= x; i /= x; return *this;}
	Complex &operator/=(const Complex &z)
		{
			double a, tmp=r;
			a = z.r*z.r + z.i*z.i;
			r = (r*z.r + i*z.i)/a;
			i = (i*z.r - tmp*z.i)/a;
			return *this;
		}



	/* conjugation, polar angle, and absolute value */
	
	friend Complex conjg(const Complex &);
	friend double arg(const Complex &);
	friend double fabs(const Complex &);
			
		

	/* arithmetic with two complex numbers */
	
	friend Complex operator+(const Complex &, const Complex &);
	friend Complex operator-(const Complex &, const Complex &);
	friend Complex operator*(const Complex &, const Complex &);
	friend Complex operator/(const Complex &, const Complex &);



	/* arithmetic with real numbers */
	
	friend Complex operator+(const Complex &, const double);
	friend Complex operator+(const double, const Complex &);
	friend Complex operator-(const Complex &, const double);
	friend Complex operator-(const double, const Complex &);
	friend Complex operator*(const Complex &, const double);
	friend Complex operator*(const double, const Complex &);
	friend Complex operator/(const double, const Complex &);
	friend Complex operator/(const Complex &, const double);



	/* exponential, logarithm, and powers */

	
	friend Complex exp(const Complex &);
	friend Complex log(const Complex &);
	friend Complex pow(const Complex &, Complex &);
	friend Complex pow(const Complex &, const double);
	friend Complex pow(const double, const Complex &);
	friend Complex sqrt(const Complex &);


	/* hyperbolic functions */
	
	friend Complex cosh(const Complex &);
	friend Complex sinh(const Complex &);
	friend Complex tanh(const Complex &);


	/* trigonometric functions */
	
	friend Complex cos(const Complex &);
	friend Complex sin(const Complex &);
	friend Complex tan(const Complex &);


	/* inverse hyperbolic functions */
	
	friend Complex acosh(const Complex &);
	friend Complex asinh(const Complex &);
	friend Complex atanh(const Complex &);


	/* inverse trigonometric functions */
	
	friend Complex acos(const Complex &);
	friend Complex asin(const Complex &);
	friend Complex atan(const Complex &);
};	
			
// inline double Re(Complex &z) {return z.r;}
// inline double Im(Complex &z) {return z.i;}
inline double Re(Complex z) {return z.r;}
inline double Im(Complex z) {return z.i;}
istream&  operator >> (istream& s, Complex& x);
ostream&  operator << (ostream& s, Complex& x);

#include <cVector.h>
		
	
#endif
