void splint1d(double *xa, double *ya, double *y2a, int n, double x, double *y, double *yd, double *ydd)
{
  int klo,khi;
  double h,b,a;
  void nrerror(char *error_text);

  klo=(int)( (x-xa[1])/(xa[n]-xa[1])*(double)(n-1) ) + 1;
  klo=klo<n ? klo : n-1;
  khi=klo+1;

  h=xa[khi]-xa[klo];
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;

  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  *yd=(-ya[klo]+ya[khi])/h + (-(3.0*a*a-1.0)*y2a[klo]+(3.0*b*b-1.0)*y2a[khi])
    *h/6.0;
  *ydd = a*y2a[klo]+b*y2a[khi];
}
