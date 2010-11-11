#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <kevin_complex.h>

using namespace std;

#ifdef CDEBUG
int Complex::nlive = 0;
#endif

KComplex conjg(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = z.r;
  tmp.i = -z.i;
  return tmp;
}

double arg(const KComplex &z)
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

double fabs(const KComplex &z)
{
  return sqrt(z.r*z.r + z.i*z.i);
}



KComplex operator+(const KComplex &z1, const KComplex &z2)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif 
  tmp.r=z1.r+z2.r; 
  tmp.i=z1.i+z2.i; 
  return tmp;
}

KComplex operator-(const KComplex &z1, const KComplex &z2)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif 
  tmp.r=z1.r-z2.r; 
  tmp.i=z1.i-z2.i; 
  return tmp;
}

KComplex operator*(const KComplex &z1, const KComplex &z2)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif 
  tmp.r=z1.r*z2.r-z1.i*z2.i;
  tmp.i = z1.r * z2.i + z1.i * z2.r;
  return tmp;
}

KComplex operator/(const KComplex &z1, const KComplex &z2)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  double a;
  a = z2.r*z2.r + z2.i*z2.i;
  tmp.r = (z1.r*z2.r + z1.i*z2.i)/a;
  tmp.i = (z1.i*z2.r - z1.r*z2.i)/a;
  return tmp;
}

KComplex operator+(const KComplex &z, const double x)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = z.r + x;
  tmp.i = z.i;
  return tmp;
}


KComplex operator+(const double x, const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = z.r + x;
  tmp.i = z.i;
  return tmp;
}
KComplex operator-(const KComplex &z, const double x)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = z.r - x;
  tmp.i = z.i;
  return tmp;
}


KComplex operator-(const double x, const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = x - z.r;
  tmp.i = -z.i;
  return tmp;
}

KComplex operator*(const KComplex &z, const double x)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = z.r*x;
  tmp.i = z.i*x;
  return tmp;
}
KComplex operator*(const double x, const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = z.r*x;
  tmp.i = z.i*x;
  return tmp;
}

KComplex operator/(const double x, const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  double a;
  a = z.r*z.r + z.i*z.i;
  tmp.r = x*z.r/a;
  tmp.i = -x*z.i/a;
  return tmp;
}

KComplex operator/(const KComplex &z, const double x)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = z.r/x;
  tmp.i = z.i/x;
  return tmp;
}

KComplex exp(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  double e = exp(z.r);
  tmp.r = e * cos(z.i);
  tmp.i = e * sin(z.i);
  return tmp;
}
KComplex log(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  double a=z.r*z.r + z.i*z.i;
  
  if (a==0.0) 
    {
      puts("domain error in KComplex log");
      exit(0);
    }
  tmp.r = 0.5*log(a);
  tmp.i = arg(z);
  return tmp;
}
KComplex pow(const KComplex &base, const KComplex &e)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp = exp(e*log(base));
  return tmp;
}
KComplex pow(const KComplex &base, const double e)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp = exp(e*log(base));
  return tmp;
}
KComplex pow(const double base, const KComplex &e)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp = exp(e*log(base));
  return tmp;
}
KComplex sqrt(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  
  if (fabs(z)==0.0) {tmp.r=0.0; tmp.i=0.0;}
  else tmp = exp(0.5*log(z));
  
  return tmp;
}


KComplex cosh(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = cosh(z.r)*cos(z.i);
  tmp.i = sinh(z.r)*sin(z.i);
  return tmp;
}
KComplex sinh(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = sinh(z.r)*cos(z.i);
  tmp.i = cosh(z.r)*sin(z.i);
  return tmp;
}
KComplex tanh(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex c;
#else
  static KComplex c;
#endif
  c = cosh(z);
  if (fabs(c)==0.0) 
    {
      puts("domain error in tanh");
      exit(0);
    }
  return sinh(z)/c;
}

KComplex asinh(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp = log(z + sqrt(1.0 + z*z));
  return tmp;
}

KComplex acosh(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp = log(z + sqrt(z*z - 1.0));
  return tmp;
}

KComplex atanh(const KComplex &z)
{
	if (fabs(z-1.0) == 0.0 || fabs(z+1.0) == 0.0)
	{
		puts("domain error in KComplex atanh()");
		exit(0);
	}
	return 0.5*log((1.0 + z)/(1.0 - z));
}


KComplex cos(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = cos(z.r)*cosh(z.i);
  tmp.i = -sin(z.r)*sinh(z.i);
  return tmp;
}


KComplex sin(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex tmp;
#else
  static KComplex tmp;
#endif
  tmp.r = sin(z.r)*cosh(z.i);
  tmp.i = cos(z.r)*sinh(z.i);
  return tmp;
}

KComplex tan(const KComplex &z)
{
#if defined _REENTRANT || defined _THREAD_SAFE
  KComplex c;
#else
  static KComplex c;
#endif
	
  c = cos(z);
  if (fabs(c)==0.0) 
    {
      puts("domain error in tan");
      exit(0);
    }
  return sin(z)/c;
}


KComplex acos(const KComplex &z)
{
  static KComplex i(0., 1.);

  return acosh(z)/i;
}

KComplex asin(const KComplex &z)
{
  static KComplex i(0., 1.);

  return asinh(i*z)/i;
}

KComplex atan(const KComplex &z)
{
  static KComplex i(0., 1.);
  
  return atanh(i*z)/i;
}


ostream& operator << (ostream& s, KComplex& x)
{
  return s << "(" << x.real() << ", " << x.imag() << ")" ;
}

istream& operator >> (istream& s, KComplex& x)
{

#if __GNUC__ && (__GNUC__ < 3)
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
  x = KComplex(r, i);
  return s;
}
