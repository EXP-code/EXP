
#ifndef _k_complex_h
#define _k_complex_h

#include <iostream>

using namespace std;

class CVector;
class CMatrix;
class KComplex;

// #define CDEBUG 1


class KComplex
{
 protected:
	double r;
	double i;
#ifdef CDEBUG
	static int nlive;
#endif
 public:
	friend class CVector;
	friend class CMatrix;

	/* constructors */
	
#ifdef CDEBUG
	KComplex(void) {r=0.0; i=0.0; nlive++;}
	KComplex(const double x) {r=x; i=0.0; nlive++;}
	KComplex(const double x, const double y) {r=x; i=y; nlive++;}
	KComplex(const KComplex &z) {r=z.r; i=z.i; nlive++;}
	~KComplex(void) {nlive--;}
#else
	KComplex(void) {r=0.0; i=0.0;}
	KComplex(const double x) {r=x; i=0.0;}
	KComplex(const double x, const double y) {r=x; i=y;}
	KComplex(const KComplex &z) {r=z.r; i=z.i;}
#endif
	

	/* access to real and imag parts */
	
	double & real(void) {return r;}
	double & imag(void) {return i;}	
	friend inline double Re(KComplex z);
	friend inline double Im(KComplex z);
	

	
	/* simple output */
	
	void print(void)	{cout << r << " + " << i << " I\n";}
#ifdef CDEBUG
	static void count(void)	{
	  cerr << "number of live complexes: " << nlive << '\n';
	}
#endif	
	
	
	/* assignment */
	
	KComplex &operator=(const KComplex &z) {r=z.r; i=z.i; return *this;}
	KComplex &operator=(const double &x) {r=x; i=0.0; return *this;}

	
	/* unary minus */
	
	KComplex operator-(void)	
	{
			KComplex tmp;
			tmp.r = -r;
			tmp.i = -i;
			return tmp;
		}
	

	
	/* reflexive arithmetic operators */
	
	KComplex &operator+=(const KComplex &z) {r+=z.r; i+=z.i; return *this;}
	KComplex &operator+=(double x) {r+=x; return *this;}
	KComplex &operator-=(const KComplex &z) {r-=z.r; i-=z.i; return *this;}
	KComplex &operator-=(double x) {r-=x; return *this;}
	KComplex &operator*=(const KComplex &z) 
	  {
			double tmp=r; 
			r=r*z.r - i*z.i; 
			i=i*z.r + tmp*z.i; 
			return *this;
		}
	KComplex &operator*=(double x) {r *= x; i *= x; return *this;}
	KComplex &operator/=(double x) {r /= x; i /= x; return *this;}
	KComplex &operator/=(const KComplex &z)
		{
		  double a, tmp=r;
			a = z.r*z.r + z.i*z.i;
			r = (r*z.r + i*z.i)/a;
			i = (i*z.r - tmp*z.i)/a;
			return *this;
		}
	

	
	/* conjugation, polar angle, and absolute value */
	
	friend KComplex conjg(const KComplex &);
	friend double arg(const KComplex &);
	friend double fabs(const KComplex &);
	
		
	
	/* arithmetic with two complex numbers */
	
	friend KComplex operator+(const KComplex &, const KComplex &);
	friend KComplex operator-(const KComplex &, const KComplex &);
	friend KComplex operator*(const KComplex &, const KComplex &);
	friend KComplex operator/(const KComplex &, const KComplex &);

	

	/* arithmetic with real numbers */
	
	friend KComplex operator+(const KComplex &, const double);
	friend KComplex operator+(const double, const KComplex &);
	friend KComplex operator-(const KComplex &, const double);
	friend KComplex operator-(const double, const KComplex &);
	friend KComplex operator*(const KComplex &, const double);
	friend KComplex operator*(const double, const KComplex &);
	friend KComplex operator/(const double, const KComplex &);
	friend KComplex operator/(const KComplex &, const double);
	

	
	/* exponential, logarithm, and powers */
	
	
	friend KComplex exp(const KComplex &);
	friend KComplex log(const KComplex &);
	friend KComplex pow(const KComplex &, KComplex &);
	friend KComplex pow(const KComplex &, const double);
	friend KComplex pow(const double, const KComplex &);
	friend KComplex sqrt(const KComplex &);
	

	/* hyperbolic functions */
	
	friend KComplex cosh(const KComplex &);
	friend KComplex sinh(const KComplex &);
	friend KComplex tanh(const KComplex &);

	
	/* trigonometric functions */
	
	friend KComplex cos(const KComplex &);
	friend KComplex sin(const KComplex &);
	friend KComplex tan(const KComplex &);
	

	/* inverse hyperbolic functions */
	
	friend KComplex acosh(const KComplex &);
	friend KComplex asinh(const KComplex &);
	friend KComplex atanh(const KComplex &);

	
	/* inverse trigonometric functions */
	
	friend KComplex acos(const KComplex &);
	friend KComplex asin(const KComplex &);
	friend KComplex atan(const KComplex &);
};	
			
// inline double Re(KComplex &z) {return z.r;}
// inline double Im(KComplex &z) {return z.i;}
inline double Re(KComplex z) {return z.r;}
inline double Im(KComplex z) {return z.i;}
istream&  operator >> (istream& s, KComplex& x);
ostream&  operator << (ostream& s, KComplex& x);

#include <CVector.h>

	
#endif
