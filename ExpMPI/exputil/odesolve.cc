#include <stdio.h>
#include <math.h>
#include <numerical.h>


void integrate_ode(
	double *x,
	double t0, double t1,
	double *dt, double eps,
	int n,
	ode_derivs derivs,
	ode_integrator integrator)
{
	double t;

	t=t0;
	
	do
	{
		*dt = MIN(*dt, t1-t);	
		onestep(x, &t, n, dt, eps, derivs, integrator);
	}
	while (t < t1);
}
		
		



void onestep(
	double *x,
	double *t,
	int n,
	double *dt, 
	double eps,
	ode_derivs derivs,
	ode_integrator integrator)
{
	double *dxdt, *xscale, hlast, hnext;
	int i;

	xscale = nr_vector(1, n);
	dxdt = nr_vector(1, n);

	derivs(*t, x, dxdt);
	for (i=1; i<=n; i++)
		xscale[i] = fabs(x[i]) + fabs(*dt * dxdt[i]) + TINY;

	integrator(x, dxdt, n, t, *dt, eps, xscale, &hlast, &hnext, derivs);

	*dt = hnext;

	free_nr_vector(xscale, 1, n);
	free_nr_vector(dxdt, 1, n);
}


void constant_step(
	double *x,
	double *t,
	int n,
	double dt, 
	ode_derivs derivs,
	ode_rk algorithm)
{
	double *dxdt;

	dxdt = nr_vector(1, n);

	derivs(*t, x, dxdt);
	algorithm(x, dxdt, n, *t, dt, x, derivs);

	*t += dt;
	free_nr_vector(dxdt, 1, n);
}






