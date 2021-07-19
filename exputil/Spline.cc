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

#include <cmath>
#include <Eigen/Eigen>

void Spline(const Eigen::VectorXd &x,
	    const Eigen::VectorXd &y,
	    double yp1, double ypn, Eigen::VectorXd &y2)
{
  double d1,d2,p,qn,un;
  double sig;
  
  int sz = x.size();
  Eigen::VectorXd u(sz-1);

/*     Boundary conditions obtained by fixing third derivative as computed
       by divided differences */
  if (yp1 < -0.99e30) {
    y2[0]=1.0;
    d2 = ((y[3]-y[2])/(x[3]-x[2]) - (y[2]-y[1])/(x[2]-x[1]))/(x[3]-x[1]);
    d1 = ((y[2]-y[1])/(x[2]-x[1]) - (y[1]-y[0])/(x[1]-x[0]))/(x[2]-x[0]);
    u[0] = -6.0*(d2-d1)*(x[1]-x[0])/(x[3]-x[0]);
  }
/*     "Normal" zero second derivative boundary conditions */
  else if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
/*      Known first derivative */
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (int i=1; i<sz-1; i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

/*     Boundary conditions obtained by fixing third derivative as computed
       by divided differences */
  if (ypn < -0.99e30) {
    d2 = ((y[sz-1]-y[sz-2])/(x[sz-1]-x[sz-2]) - 
	  (y[sz-2]-y[sz-3])/(x[sz-2]-x[sz-3]))/(x[sz-1]-x[sz-3]);
    d1 = ((y[sz-2]-y[sz-3])/(x[sz-2]-x[sz-3]) - 
	  (y[sz-3]-y[sz-4])/(x[sz-3]-x[sz-4]))/(x[sz-2]-x[sz-4]);
    qn = -1.0;
    un = 6.0*(d2-d1)*(x[sz-1]-x[sz-2])/(x[sz-1]-x[sz-4]);
  }
/*     "Normal" zero second derivative boundary conditions */
  else if (ypn > 0.99e30)
    qn=un=0.0;
/*      Known first derivative */
  else {
    qn=0.5;
    un=(3.0/(x[sz-1]-x[sz-2]))*(ypn-(y[sz-1]-y[sz-2])/(x[sz-1]-x[sz-2]));
  }
  y2[sz-1]=(un-qn*u[sz-2])/(qn*y2[sz-2]+1.0);
  for (int k=sz-2; k>=0; k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
}

