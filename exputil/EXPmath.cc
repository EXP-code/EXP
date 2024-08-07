#include <iostream>
#include <iomanip>
#include <cmath>

#include "EXPmath.H"
#include "laguerre_polynomial.hpp"

namespace AltMath
{

  double bessj1(double x)
  {
    double ax,z;
    double xx,y,ans,ans1,ans2;
  
    if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
						+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
      ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
					     +y*(99447.43394+y*(376.9991397+y*1.0))));
      ans=ans1/ans2;
    } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
				 +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
			    +y*(0.8449199096e-5+y*(-0.88228987e-6
						   +y*0.105787412e-6)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
      if (x < 0.0) ans = -ans;
    }
    return ans;
  }
  
  double bessj0(double x)
  {
    double ax,z;
    double xx,y,ans,ans1,ans2;
    
    if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
					      +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
					    +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
    } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
				      +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
				 +y*(-0.6911147651e-5+y*(0.7621095161e-6
							 -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    }
    return ans;
  }

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

  double bessj(int n,double x)
  {
    int j,jsum,m;
    double ax,bj,bjm,bjp,sum,tox,ans;
    
    if (n < 2) {
      std::cerr << "Index n less than 2 in BESSJ" << std::endl;
      exit(-1);
    }
  
    ax=fabs(x);
    if (ax == 0.0)
      return 0.0;
    else if (ax > (double) n) {
      tox=2.0/ax;
      bjm=bessj0(ax);
      bj=bessj1(ax);
      for (j=1;j<n;j++) {
	bjp=j*tox*bj-bjm;
	bjm=bj;
	bj=bjp;
      }
      ans=bj;
    } else {
      tox=2.0/ax;
      m=2*((n+(int) sqrt(ACC*n))/2);
      jsum=0;
      bjp=ans=sum=0.0;
      bj=1.0;
      for (j=m;j>0;j--) {
	bjm=j*tox*bj-bjp;
	bjp=bj;
	bj=bjm;
	if (fabs(bj) > BIGNO) {
	  bj *= BIGNI;
	  bjp *= BIGNI;
	  ans *= BIGNI;
	  sum *= BIGNI;
	}
	if (jsum) sum += bj;
	jsum=!jsum;
	if (j == n) ans=bjp;
      }
      sum=2.0*sum-bj;
      ans /= sum;
    }
    return  x < 0.0 && n%2 == 1 ? -ans : ans;
  }
  
#undef ACC
#undef BIGNO
#undef BIGNI

  double bessk0(double x)
  {
    double y,ans;
    double bessi0(double);

    if (x <= 2.0) {
      y=x*x/4.0;
      ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
						  +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
										    +y*(0.10750e-3+y*0.74e-5))))));
    } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
					   +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
										+y*(-0.251540e-2+y*0.53208e-3))))));
    }
    return ans;
  }

  double bessk(int n,double x)
  {
    int j;
    double bk,bkm,bkp,tox;
    double bessk0(double),bessk1(double);
    
    if (n < 2) throw std::runtime_error("Index n less than 2 in BESSK");
    tox=2.0/x;
    bkm=bessk0(x);
    bk=bessk1(x);
    for (j=1;j<n;j++) {
      bkp=bkm+j*tox*bk;
      bkm=bk;
      bk=bkp;
    }
    return bk;
  }

  double bessk1(double x)
  {
    double y,ans;
    double bessi1(double);
    
    if (x <= 2.0) {
      y=x*x/4.0;
      ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
						 +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
										   +y*(-0.110404e-2+y*(-0.4686e-4)))))));
    } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
					   +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
										+y*(0.325614e-2+y*(-0.68245e-3)))))));
    }
    return ans;
  }

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

  double bessi(int n,double x)
  {
    int j;
    double bi,bim,bip,tox,ans;
    double bessi0(double );
    
    if (n < 2) std::runtime_error("Index n less than 2 in BESSI");
    if (x == 0.0)
      return 0.0;
    else {
      tox=2.0/fabs(x);
      bip=ans=0.0;
      bi=1.0;
      for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
	bim=bip+j*tox*bi;
	bip=bi;
	bi=bim;
	if (fabs(bi) > BIGNO) {
	  ans *= BIGNI;
	  bi *= BIGNI;
	  bip *= BIGNI;
	}
	if (j == n) ans=bip;
      }
      ans *= bessi0(x)/bi;
      return  x < 0.0 && n%2 == 1 ? -ans : ans;
    }
  }

#undef ACC
#undef BIGNO
#undef BIGNI

  double bessi0(double x)
  {
    double ax,ans;
    double y;
    
    if ((ax=fabs(x)) < 3.75) {
      y=x/3.75;
      y*=y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
					   +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
					    +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
									       +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
														    +y*0.392377e-2))))))));
    }
    return ans;
  }

  
  double bessi1(double x)
  {
    double ax,ans;
    double y;
    
    if ((ax=fabs(x)) < 3.75) {
      y=x/3.75;
      y*=y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
						 +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
    } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
					   -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
					 +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
  }

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

  double jn_sph(int n, double x)
  {
    int j,m;
    double ax,bj,bjm,bjp,tox,ans;
    double jn_sph0(double x),jn_sph1(double x),jn_sphm1(double x);
    
    switch (n) {
    case -1:
      return (jn_sphm1(x));
    case 0:
      return (jn_sph0(x));
    case 1:
      return (jn_sph1(x));
    }
    
    ax=fabs(x);
    if (ax == 0.0)
      return 0.0;
    else if (ax > 0.5+ (double) n) {
      tox=1.0/ax;
      bjm=jn_sph0(ax);
      bj=jn_sph1(ax);
      for (j=1;j<n;j++) {
	bjp=(2*j+1)*tox*bj-bjm;
	bjm=bj;
	bj=bjp;
      }
      ans=bj;
    } else {
      tox=1.0/ax;
      m=2*((n+(int) sqrt(ACC*n))/2);
      bjp=ans=0.0;
      bj=1.0;
      for (j=m;j>0;j--) {
	bjm=(2*j+1)*tox*bj-bjp;
	bjp=bj;
	bj=bjm;
	if (fabs(bj) > BIGNO) {
	  bj *= BIGNI;
	  bjp *= BIGNI;
	  ans *= BIGNI;
	}
	if (j == n) ans=bjp;
      }
      ans *= jn_sph0(ax)/bj;
    }
    return  x < 0.0 && n%2 == 1 ? -ans : ans;
  }
  
#undef ACC
#undef BIGNO
#undef BIGNI

#define TOL 1.0e-7
  double jn_sph0(double x)
  {
    if (fabs(x)<TOL)
      return 1.0-x*x/6.0;
    else
      return sin(x)/x;
  }

  double jn_sph1(double x)
  {
    if (fabs(x)<TOL)
      return x/3.0;
    else
      return (sin(x)/x - cos(x))/x;
  }

  double jn_sphm1(double x)
  {
    return cos(x)/x;
  }
#undef TOL

  double check_nu(double nu)
  {
    int n = floor(nu);
    double nu1 = n;
    if (fabs(n - nu1) > 1.0e-16)
      throw std::runtime_error("AltMath: Bessel implementation is only valid for integer-order arguments");
    return n;
  }

  double sph_bessel(unsigned n, double x) { return jn_sph(n, x); }

  double cyl_bessel_j(double nu, double x)
  {
    double bessj(int, double), bessj0(double), bessj1(double);
  
    int n = check_nu(nu);
    if (n==0) return bessj0(x);
    if (n==1) return bessj1(x);
    return bessj(n, x);
  }

  double cyl_bessel_k(double nu, double x)
  {
    double bessk(int, double), bessk0(double), bessk1(double);

    int n = check_nu(nu);
    if (n==0) return bessk0(x);
    if (n==1) return bessk1(x);
    return bessk(n, x);
  }

  double cyl_bessel_i(double nu, double x)
  {
    double bessi(int, double), bessi0(double), bessi1(double);

    int n = check_nu(nu);
    if (n==0) return bessi0(x);
    if (n==1) return bessi1(x);
    return bessi(n, x);
  }

  double assoc_laguerre(unsigned l, unsigned n, double x)
  {
    double xi[] = {x};
    auto yi = l_polynomial ((int)l, (int)n, xi);
    double ret = yi[0]; delete [] yi;
    return ret;
  }
  
}
