#include <math.h>

#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);

int frprmn(double *p, int n, double ftol, int *iter, double *fret, double (*func) (/* ??? */), void (*dfunc) (/* ??? */))
{
	int j,its,ITMAX;
	double gg,gam,fp,dgg;
	double *g,*h,*xi,*dvector();
	void linmin(double *p, double *xi, int n, double *fret, double (*func) (/* ??? */)),nrerror(),free_dvector();

	ITMAX = *iter;

	g=dvector(1,n);
	h=dvector(1,n);
	xi=dvector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return 0;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
/*		  dgg += xi[j]*xi[j];	*/
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return 0;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	return its;
}

#undef EPS
#undef FREEALL
