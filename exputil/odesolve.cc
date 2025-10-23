#include <iostream>
#include <iomanip>
#include <cmath>
#include "numerical.H"

void integrate_ode(Eigen::VectorXd& x,
		   double t0, double t1, double &dt, double eps,
		   int n, ode_derivs derivs, ode_integrator integrator)
{
  double t = t0;
	
  do
    {
      dt = std::min<double>(dt, t1-t);	
      onestep(x, t, n, dt, eps, derivs, integrator);
    }
  while (t < t1);
}



void onestep(
	     Eigen::VectorXd& x,
	     double& t,
	     int n,
	     double& dt, 
	     double eps,
	     ode_derivs derivs,
	     ode_integrator integrator)
{
  Eigen::VectorXd xscale(n);
  Eigen::VectorXd dxdt(n);
  double hlast, hnext;

  derivs(t, x, dxdt);
  for (int i=0; i<n; i++)
    xscale[i] = fabs(x[i]) + fabs(dt * dxdt[i]) + TINY;

  integrator(x, dxdt, n, t, dt, eps, xscale, hlast, hnext, derivs);

  dt = hnext;
}


void constant_step(
		   Eigen::VectorXd& x,
		   double& t,
		   int n,
		   double dt, 
		   ode_derivs derivs,
		   ode_rk algorithm)
{
    Eigen::VectorXd dxdt(n);

    derivs(t, x, dxdt);
    algorithm(x, dxdt, n, t, dt, x, derivs);

    t += dt;
}
