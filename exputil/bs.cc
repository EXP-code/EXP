#include <iostream>
#include <iomanip>
#include <cmath>
#include "numerical.H"

#define IMAX 11
#define NUSE 7
#define SHRINK 0.95
#define GROW 1.2

static Eigen::MatrixXd d;	/* defining declaration */
static Eigen::VectorXd x;

void bsstep(
	    Eigen::VectorXd& y,
	    Eigen::VectorXd& dydx,
	    int nv,
	    double& xx,
	    double htry,
	    double eps,
	    Eigen::VectorXd& yscal,
	    double& hdid,
	    double& hnext,
	    ode_derivs derivs)
{
  static int nseq[IMAX+1]={0,2,4,6,8,12,16,24,32,48,64,96};
  void mmid(Eigen::VectorXd&, Eigen::VectorXd&, int, double, double, int,
	    Eigen::VectorXd&, ode_derivs);
  void rzextr(int, double, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&, int, int);

  Eigen::VectorXd ysav(nv);
  Eigen::VectorXd dysav(nv);
  Eigen::VectorXd yseq(nv);
  Eigen::VectorXd yerr(nv);
  Eigen::VectorXd x(IMAX);
  Eigen::MatrixXd d(nv, NUSE);

  double h=htry;
  double xsav=xx;
  for (int i=0;i<nv;i++) {
    ysav[i]=y[i];
    dysav[i]=dydx[i];
  }

  for (;;) {
    for (int i=0; i<IMAX; i++) {
      mmid(ysav, dysav, nv, xsav, h, nseq[i], yseq, derivs);
      double temp=h/nseq[i];
      double xest=temp*temp;
      rzextr(i,xest, yseq, y, yerr, nv, NUSE);

      double errmax=0.0;
      for (int j=0;j<nv;j++)
	if (errmax < fabs(yerr[j]/yscal[j]))
	  errmax=fabs(yerr[j]/yscal[j]);
      errmax /= eps;
      if (errmax < 1.0) {
	xx += h;
	hdid=h;
	hnext = i==NUSE? h*SHRINK : i==NUSE-1?
	  h*GROW : (h*nseq[NUSE-1])/nseq[i];
	return;
      }
    }
    h *= 0.25;
    for (int i=0; i<(IMAX-NUSE)/2;i++) h /= 2.0;
    if ((xx+h) == (xx)) {
      std::cerr << "Step size underflow in BSSTEP" << std::endl;
      exit(-1);
    }
  }
}

#undef IMAX
#undef NUSE
#undef SHRINK
#undef GROW

void mmid(
	  Eigen::VectorXd& y,
	  Eigen::VectorXd& dydx,
	  int nvar,
	  double xs,
	  double htot,
	  int nstep,
	  Eigen::VectorXd& yout,
	  ode_derivs derivs)
{
  Eigen::VectorXd ym(nvar);
  Eigen::VectorXd yn(nvar);

  double h = htot/nstep;
  for (int i=0; i<nvar; i++) {
    ym[i]=y[i];
    yn[i]=y[i]+h*dydx[i];
  }

  double x=xs+h;
  derivs(x,yn,yout);

  double h2=2.0*h;
  for (int n=2; n<=nstep; n++) {
    for (int i=0;i<nvar;i++) {
      double swap = ym[i] + h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;

    derivs(x,yn,yout);
  }
  for (int i=0; i<nvar;i++)
    yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
}



void rzextr(
	int iest,
	double xest,
	Eigen::VectorXd& yest,
	Eigen::VectorXd& yz,
	Eigen::VectorXd& dy,
	int nv,
	int nuse)
{
  Eigen::VectorXd fx(nuse);

  x[iest] = xest;
  if (iest == 0)
    for (int j=0; j<nv;j++) {
      yz[j]   = yest[j];
      d(j, 0) = yest[j];
      dy[j]   = yest[j];
    }
  else {
    int m1=(iest < nuse ? iest : nuse);
    for (int k=0;k<m1-1; k++)
      fx[k+1]=x[iest-k]/xest;

    for (int j=0;j<nv;j++) {
      double yy = yest[j], ddy;
      double v  = d(j, 0);
      double c  = yy;
      d(j, 0)   = yy;
      for (int k=0; k<m1; k++) {
	double b1 = fx[k]*v;
	double b  = b1-c;
	if (b) {
	  b = (c-v)/b;
	  ddy = c*b;
	  c = b1*b;
	} else
	  ddy = v;
	if (k != m1) v = d(j, k);
	d(j, k) = ddy;
	yy += ddy;
      }
      dy[j] = ddy;
      yz[j] = yy;
    }
  }
}


