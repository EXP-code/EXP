#include <stdio.h>
#include <math.h>
#include <numerical.h>


#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666		/* 1.0/15.0 */
#define SAFETY 0.9
#define ERRCON 6.0e-4


void rkqc(	
	double *y, 
	double *dydx, 
	int n, 
	double *x,
	double htry, 
	double eps,
	double *yscal,
	double *hdid,
	double *hnext,
	ode_derivs derivs)
{
	int i;
	double xsav,hh,h,temp,errmax;
	double *dysav,*ysav,*ytemp;

	dysav=nr_vector(1,n);
	ysav=nr_vector(1,n);
	ytemp=nr_vector(1,n);
	xsav=(*x);
	for (i=1;i<=n;i++) {
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	h=htry;
	for (;;) {
		hh=0.5*h;
		rk4(ysav,dysav,n,xsav,hh,ytemp,derivs);
		*x=xsav+hh;
		(*derivs)(*x,ytemp,dydx);
		rk4(ytemp,dydx,n,*x,hh,y,derivs);
		*x=xsav+h;
		if (*x == xsav) nrerror("Step size too small in routine RKQC");
		rk4(ysav,dysav,n,xsav,h,ytemp,derivs);
		errmax=0.0;
		for (i=1;i<=n;i++) {
			ytemp[i]=y[i]-ytemp[i];
			temp=fabs(ytemp[i]/yscal[i]);
			if (errmax < temp) errmax=temp;
		}
		errmax /= eps;
		if (errmax <= 1.0) {
			*hdid=h;
			*hnext=(errmax > ERRCON ?
				SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
			break;
		}
		h=SAFETY*h*exp(PSHRNK*log(errmax));
	}
	for (i=1;i<=n;i++) y[i] += ytemp[i]*FCOR;
	free_nr_vector(ytemp,1,n);
	free_nr_vector(dysav,1,n);
	free_nr_vector(ysav,1,n);
}

#undef PGROW
#undef PSHRNK
#undef FCOR
#undef SAFETY
#undef ERRCON



void rk4(
	double *y,
	double *dydx,
	int n,
	double x,
	double h,
	double *yout,
	ode_derivs derivs)
{
	int i;
	double xh,hh,h6,*dym,*dyt,*yt;

	dym=nr_vector(1,n);
	dyt=nr_vector(1,n);
	yt=nr_vector(1,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=1;i<=n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_nr_vector(yt,1,n);
	free_nr_vector(dyt,1,n);
	free_nr_vector(dym,1,n);
}

