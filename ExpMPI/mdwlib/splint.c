void splint(double *xa, double *ya, double *y2a, int n, double x, double *y, double *yd, double *ydd)
{
	int klo,khi,k;
	double h,b,a;
	void nrerror(char *);

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
   *yd=(-ya[klo]+ya[khi])/h +
		(-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
     	*h/6.0;
	*ydd = a*y2a[klo]+b*y2a[khi];
}
