void splint1(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
  int klo,khi;
  double h,b,a;
  void nrerror();

  klo=(int)( (x-xa[1])/(xa[n]-xa[1])*(double)(n-1) ) + 1;
  klo=klo<n ? klo : n-1;
  khi=klo+1;
  h=xa[khi]-xa[klo];
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
