#include <math.h>
#include <numerical.h>

#define IMAX 11
#define NUSE 7
#define SHRINK 0.95
#define GROW 1.2

static double **d=0,*x=0;	/* defining declaration */

void bsstep(
	double *y,
	double *dydx,
	int nv,
	double *xx,
	double htry,
	double eps,
	double *yscal,
	double *hdid,
	double *hnext,
	ode_derivs derivs)
{
	int i,j;
	double xsav,xest,h,errmax,temp;
	double *ysav,*dysav,*yseq,*yerr;
	static int nseq[IMAX+1]={0,2,4,6,8,12,16,24,32,48,64,96};
	void mmid(double *, double *, int, double, double, int, double *,
		ode_derivs);
	void rzextr(int, double, double *, double *, double *, int, int);

	ysav=nr_vector(1,nv);
	dysav=nr_vector(1,nv);
	yseq=nr_vector(1,nv);
	yerr=nr_vector(1,nv);
	x=nr_vector(1,IMAX);
	d=nr_matrix(1,nv,1,NUSE);
	h=htry;
	xsav=(*xx);
	for (i=1;i<=nv;i++) {
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	for (;;) {
		for (i=1;i<=IMAX;i++) {
			mmid(ysav,dysav,nv,xsav,h,nseq[i],yseq,derivs);
			xest=(temp=h/nseq[i],temp*temp);
			rzextr(i,xest,yseq,y,yerr,nv,NUSE);
			errmax=0.0;
			for (j=1;j<=nv;j++)
				if (errmax < fabs(yerr[j]/yscal[j]))
					errmax=fabs(yerr[j]/yscal[j]);
			errmax /= eps;
			if (errmax < 1.0) {
				*xx += h;
				*hdid=h;
				*hnext = i==NUSE? h*SHRINK : i==NUSE-1?
					h*GROW : (h*nseq[NUSE-1])/nseq[i];
				free_nr_matrix(d,1,nv,1,NUSE);
				free_nr_vector(x,1,IMAX);
				free_nr_vector(yerr,1,nv);
				free_nr_vector(yseq,1,nv);
				free_nr_vector(dysav,1,nv);
				free_nr_vector(ysav,1,nv);
				return;
			}
		}
		h *= 0.25;
		for (i=1;i<=(IMAX-NUSE)/2;i++) h /= 2.0;
		if ((*xx+h) == (*xx)) nrerror("Step size underflow in BSSTEP");
	}
}

#undef IMAX
#undef NUSE
#undef SHRINK
#undef GROW

void mmid(
	double *y,
	double *dydx,
	int nvar,
	double xs,
	double htot,
	int nstep,
	double *yout,
	ode_derivs derivs)
{
	int n,i;
	double x,swap,h2,h,*ym,*yn;

	ym=nr_vector(1,nvar);
	yn=nr_vector(1,nvar);
	h=htot/nstep;
	for (i=1;i<=nvar;i++) {
		ym[i]=y[i];
		yn[i]=y[i]+h*dydx[i];
	}
	x=xs+h;
	(*derivs)(x,yn,yout);
	h2=2.0*h;
	for (n=2;n<=nstep;n++) {
		for (i=1;i<=nvar;i++) {
			swap=ym[i]+h2*yout[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		x += h;
		(*derivs)(x,yn,yout);
	}
	for (i=1;i<=nvar;i++)
		yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
	free_nr_vector(yn,1,nvar);
	free_nr_vector(ym,1,nvar);
}



void rzextr(
	int iest,
	double xest,
	double *yest,
	double *yz,
	double *dy,
	int nv,
	int nuse)
{
	int m1,k,j;
	double yy,v,ddy,c,b1,b,*fx;

	fx=nr_vector(1,nuse);
	x[iest]=xest;
	if (iest == 1)
		for (j=1;j<=nv;j++) {
			yz[j]=yest[j];
			d[j][1]=yest[j];
			dy[j]=yest[j];
		}
	else {
		m1=(iest < nuse ? iest : nuse);
		for (k=1;k<=m1-1;k++)
			fx[k+1]=x[iest-k]/xest;
		for (j=1;j<=nv;j++) {
			yy=yest[j];
			v=d[j][1];
			c=yy;
			d[j][1]=yy;
			for (k=2;k<=m1;k++) {
				b1=fx[k]*v;
				b=b1-c;
				if (b) {
					b=(c-v)/b;
					ddy=c*b;
					c=b1*b;
				} else
					ddy=v;
				if (k != m1) v=d[j][k];
				d[j][k]=ddy;
				yy += ddy;
			}
			dy[j]=ddy;
			yz[j]=yy;
		}
	}
	free_nr_vector(fx,1,nuse);
}


