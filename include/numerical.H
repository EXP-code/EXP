#ifndef _numerical_H
#define _numerical_H

#include <Eigen/Eigen>

// some simple macros (some very old stuff, not yet replaced)

#ifndef MAX
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#endif
#ifndef MIN
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#endif

#define FLOAT(X)  ((float) (X))
#define DOUBLE(X) ((double) (X))
#define SQR(X) ((X)*(X))
#ifndef TINY
#define TINY 1.e-10
#endif
#define Sign(X) (((X)>0) ? (1.0) : (-1.0))

// special functions

// double BesselJ(int, double);
// double BesselI(int, double);
// double BesselK(int, double);
// double gammln(double);


const double Pi = 3.141592653589793;


double poly(double, int, double *);

// typedefs for 1D functions used in quadrature, minimization, and roots

typedef double (*func_1d)  (double);
typedef double (*func_2d)  (double, double);
typedef void   (*dfunc_1d) (double, double&, double&);


// numerical quadrature
//
void   gauleg(double, double, Eigen::VectorXd&, Eigen::VectorXd&, int);
double qadapt(double, double, func_1d, double);
double qadapt2d(double, double, double, double, func_2d, double);



// root finding
//
double zbrent(func_1d, double, double, double);

double rtbis(func_1d, double, double, double);

double rtsafe(dfunc_1d, double, double, double);

int zbrac(func_1d, Eigen::VectorXd&, Eigen::VectorXd&);




// 1D minimization
//
double brent(double, double, double, func_1d, double, double&);
void mnbrak(Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&, 
	Eigen::VectorXd&, func_1d);



// ODE solvers
//
// first, some typedefs for passing function pointers
//
typedef void (*ode_derivs)(double, Eigen::VectorXd&, Eigen::VectorXd&);

typedef void (*ode_integrator)(Eigen::VectorXd&, Eigen::VectorXd&, int,
			       double&, double, double, Eigen::VectorXd&,
			       double&, double&, ode_derivs);

typedef void (*ode_rk)(Eigen::VectorXd&, Eigen::VectorXd&, int, double, double, 
		       Eigen::VectorXd&, ode_derivs);


typedef void (*symp_derivs)(double, Eigen::VectorXd&, Eigen::VectorXd&,
			    Eigen::VectorXd&);

typedef void (*symp_integrator)(Eigen::VectorXd&, Eigen::VectorXd&,
				Eigen::VectorXd&, Eigen::VectorXd&, 
				double, double, symp_derivs);


// prototypes for ODE solvers
//
void integrate_ode(Eigen::VectorXd&, double, double,
		   double&, double, int,
		   ode_derivs, ode_integrator);

void onestep(Eigen::VectorXd&, double&, int, double&,
	     double, ode_derivs, ode_integrator);
	
void rkqc(Eigen::VectorXd&, Eigen::VectorXd&, int, double&, double, double, 
	  Eigen::VectorXd&, double&, double&, ode_derivs);
		
void rk2qc (Eigen::VectorXd&, Eigen::VectorXd&, int, double&, double, double, 
	    Eigen::VectorXd&, double&, double&, ode_derivs);

void bsstep (Eigen::VectorXd&, Eigen::VectorXd&, int, double&, 
	     double, double, Eigen::VectorXd&, double&, double&, ode_derivs);

void constant_step (Eigen::VectorXd&, Eigen::VectorXd&,
		    int, double, ode_derivs, ode_rk);

void rk4(Eigen::VectorXd&, Eigen::VectorXd&, int, double, double,
	 Eigen::VectorXd&, ode_derivs);

void sia4(Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&,
	  Eigen::VectorXd&, double, double, symp_derivs);


// interpolation stuff
//
void locate (Eigen::VectorXd&, int, double, int *);

#endif
