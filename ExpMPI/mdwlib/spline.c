/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine finds the cubic spline coefficients.  The boundary condi-
 *  tions may be one of the following three:
 *       1) "natural," that is, zero second derivatives
 *       2) first derivatives specified
 *       3) third derivatives computed from supplied data
 *
 *
 *  Call sequence:
 *  -------------
 *  void spline(x,y,n,yp1,ypn,y2);
 *
 *  int n;
 *  double x[n],y[n],yp1,ypn,y2[n];
 *
 *  Parameters:
 *  ----------
 *
 *  n        number of supplied grid points
 *  x        abcissa array
 *  y        ordinate array
 *  yp1      boundary condition at j=1
 *  ypn      boundary condition at j=n
 *  y2       array to contain spline coefficients
 *
 *  Returns:
 *  -------
 *
 *  None, spline coefficients returned by pointer
 *
 *  Notes:
 *  -----
 *
 *  If     yp1,yp2 >  1.0e30  boundary conditions (1) natural splines are used
 *  If     yp1,yp2 < -1.0e30  boundary conditions (3) approx. 3rd derivs used
 *  Otherwise                 boundary conditions (2) explicit 2nd derivs used
 *
 *  By:
 *  --
 *  Adopted from Numerical Recipes, Press et al.
 *  Third deriv. boundary condition ---  MDW 11/13/88
 *
 ***************************************************************************/

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
  int i,k;
  double d1,d2,p,qn,sig,un,*u,*dvector(int, int);
  void free_dvector(double *, int, int);
  
  u=dvector(1,n-1);

/*     Boundary conditions obtained by fixing third derivative as computed
       by divided differences */
  if (yp1 < -0.99e30) {
    y2[1]=1.0;
    d2 = ((y[4]-y[3])/(x[4]-x[3]) - (y[3]-y[2])/(x[3]-x[2]))/(x[4]-x[2]);
    d1 = ((y[3]-y[2])/(x[3]-x[2]) - (y[2]-y[1])/(x[2]-x[1]))/(x[3]-x[1]);
    u[1] = -6.0*(d2-d1)*(x[2]-x[1])/(x[4]-x[1]);
  }
/*     "Normal" zero second derivative boundary conditions */
  else if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
/*      Known first derivative */
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

/*     Boundary conditions obtained by fixing third derivative as computed
       by divided differences */
  if (ypn < -0.99e30) {
    d2 = ((y[n]-y[n-1])/(x[n]-x[n-1]) - 
	  (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))/(x[n]-x[n-2]);
    d1 = ((y[n-1]-y[n-2])/(x[n-1]-x[n-2]) - 
	  (y[n-2]-y[n-3])/(x[n-2]-x[n-3]))/(x[n-1]-x[n-3]);
    qn = -1.0;
    un = 6.0*(d2-d1)*(x[n]-x[n-1])/(x[n]-x[n-3]);
  }
/*     "Normal" zero second derivative boundary conditions */
  else if (ypn > 0.99e30)
    qn=un=0.0;
/*      Known first derivative */
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_dvector(u,1,n-1);
}

