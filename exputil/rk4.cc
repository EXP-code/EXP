#include <iostream>
#include <iomanip>
#include <cmath>

#include "numerical.H"

#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666		/* 1.0/15.0 */
#define SAFETY 0.9
#define ERRCON 6.0e-4


void rkqc(	
	  Eigen::VectorXd& y, 
	  Eigen::VectorXd& dydx, 
	  int n, 
	  double &x,
	  double htry, 
	  double eps,
	  Eigen::VectorXd& yscal,
	  double &hdid,
	  double &hnext,
	  ode_derivs derivs)
{
  Eigen::VectorXd dysav(n);
  Eigen::VectorXd ysav(n);
  Eigen::VectorXd ytemp(n);

  double xsav = x;

  for (int i=0;i<n;i++) {
    ysav[i]  = y[i];
    dysav[i] = dydx[i];
  }
  double h=htry;
  
  for (;;) {
    double hh = 0.5*h;

    rk4(ysav,dysav,n,xsav,hh,ytemp,derivs);
    x = xsav + hh;
    derivs(x,ytemp,dydx);
    rk4(ytemp,dydx,n,x,hh,y,derivs);
    x = xsav + h;
    if (x == xsav) {
      std::cerr << "Step size too small in routine RKQC" << std::endl;
      exit(-1);
    }
    rk4(ysav,dysav,n,xsav,h,ytemp,derivs);
    double errmax=0.0;
    for (int i=0; i<n;i++) {
      ytemp[i]=y[i]-ytemp[i];
      double temp=fabs(ytemp[i]/yscal[i]);
      if (errmax < temp) errmax=temp;
    }
    errmax /= eps;
    if (errmax <= 1.0) {
      hdid=h;
      hnext=(errmax > ERRCON ?
	     SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
      break;
    }
    h = SAFETY*h*exp(PSHRNK*log(errmax));
  }

  for (int i=0; i<n; i++) y[i] += ytemp[i]*FCOR;
}

#undef PGROW
#undef PSHRNK
#undef FCOR
#undef SAFETY
#undef ERRCON



void rk4(Eigen::VectorXd& y,
	 Eigen::VectorXd& dydx,
	 int n,
	 double x,
	 double h,
	 Eigen::VectorXd& yout,
	 ode_derivs derivs)
{
  Eigen::VectorXd dym(n);
  Eigen::VectorXd dyt(n);
  Eigen::VectorXd yt(n);

  double hh=h*0.5;
  double h6=h/6.0;
  double xh=x+hh;
  for (int i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
  derivs(xh,yt,dyt);

  for (int i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
  derivs(xh,yt,dym);

  for (int i=0;i<n;i++) {
    yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  derivs(x+h,yt,dyt);

  for (int i=0;i<n;i++)
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

