/**************************************************************************
*
*		Computes spherical Bessel functions (first kind) j_n(x)
*
*		From recursion relation using Miller's "downward" algorithm
*		for x<n
*
*               good for -1 <= n 
*
*		MDW: 12/26/1987
**************************************************************************/
#include <math.h>

/* test routine
main()
{
	double sbessj(),x;
	int n;

	printf("Test of j_n(x) . . .");
	printf("n,x: ");
	scanf("%d %lf",&n,&x);
	printf("j_%d(%.3le)=%.10le\n",n,x,sbessj(n,x));
}

  end of test */


#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

double sbessj(int n, double x)
{
	int j,m;
	double ax,bj,bjm,bjp,tox,ans;
	double sbessj0(double x),sbessj1(double x),sbessjm1(double x);
	void myerror();

	switch (n) {
	        case -1:
	                return (sbessjm1(x));
		case 0:
			return (sbessj0(x));
		case 1:
			return (sbessj1(x));
		}

	ax=fabs(x);
	if (ax == 0.0)
		return 0.0;
	else if (ax > 0.5+ (double) n) {
		tox=1.0/ax;
		bjm=sbessj0(ax);
		bj=sbessj1(ax);
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
		ans *= sbessj0(ax)/bj;
	}
	return  x < 0.0 && n%2 == 1 ? -ans : ans;
}

#undef ACC
#undef BIGNO
#undef BIGNI

#define TOL 1.0e-7
double sbessj0(double x)
{
	if (fabs(x)<TOL)
		return 1.0-x*x/6.0;
	else
		return sin(x)/x;
}

double sbessj1(double x)
{
	if (fabs(x)<TOL)
		return x/3.0;
	else
		return (sin(x)/x - cos(x))/x;
}

double sbessjm1(double x)
{
        return cos(x)/x;
}
#undef TOL

