#include <math.h>
#include <stdio.h>

/*
*	Test romberg routine
*/
/*
main()
{
	double	a,b,func(),rombe2(),exact();
	int		n;

	while (1)
		{
		printf("a,b,n: ");
		if (scanf("%lf %lf %d",&a,&b,&n) != 3)
			{
			printf("wrong number of arguments\n\n");
			break;
			}

		printf("ans = %35.20e   exact = %35.20e\n",rombe2(a,b,func,n),
		       exact(a,b));
		}

}

double	func(x)
double	x;
{
	return(x*x*x + x*x + x + 5.0);
}
#define  integral(z) (0.25*z*z*z*z + z*z*z/3.0 + 0.5*z*z + 5.0*z)
double	exact(a,b)
double	a,b;
{
	return(integral(b) - integral(a));
}
END OF TEST */

/*
*	Computes integral using Simpson's rule with Romberg's
*	correction with a supplied integrand routine
*
*	Synopsis:
*		answer = rombe2(a,b,f,n);
*		double	answer,a,b,f(),tol;
*		int	n;
*
* 	Use:
* 		a	initial pt.
* 		b	final pt.
* 		f	function
* 		ans	contains answer on return
* 		t	holds first column of Romberg tableux
* 
*/
#define		NROMB 12
#define         DMACH 1.0e-8

double	rombe2(double a, double b, double (*f) (double), int n)
{
	double		t[NROMB],xmesh,d,ends,x2,x4;
	int			nmax,num,i,j,k,index;

	nmax = (int)(pow(2.0,(double)n) + DMACH);
	xmesh = (b-a)/nmax;

	d = (b-a)/2.0;
	ends = (*f)(a) + (*f)(b);
	num = nmax/2;
	x2 = (*f)(a+xmesh*num);
	t[0] = d*(ends+4.0*x2)/3.0;
	for(i=2; i<=n; i++)
		{
		x4 = 0.0;
		d = d/2.0;
		index = num/2;
		for(j=1; j<=nmax/num; j++)
			{
			x4 = x4 + (*f)(a+xmesh*index);
			index = index + num;
			}
		t[i-1] = d*(ends+2.0*x2+4.0*x4)/3.0;

		x2 = x2+x4;
		num = num/2;

		}

	num=3;
	for(j=1; j <= n-1; j++)
		{
		for(k=1; k <= n-j; ++k)
			t[k-1] = t[k]+(t[k]-t[k-1])/num;
		num = 4*num+3;
		}
	return(t[0]);

}
#undef          DMACH
