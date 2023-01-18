#ifndef _numerical_H
#define _numerical_H

#include <functional>
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

// numerical quadrature
//
void   gauleg(double, double, Eigen::VectorXd&, Eigen::VectorXd&, int);
double qadapt(double, double, std::function<double(double)>, double);
double qadapt2d(double, double, double, double,
		std::function<double(double, double)>, double);



// root finding
//
double zbrent(std::function<double(double)>, double, double, double);

double rtbis(std::function<double(double)>, double, double, double);

double rtsafe(std::function<void(double, double&, double&)>, double, double, double);

int zbrac(std::function<double(double)>, Eigen::VectorXd&, Eigen::VectorXd&);




// 1D minimization
//
double brent(double, double, double, std::function<double(double)>,
	     double, double&);
void mnbrak(Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&,
	    Eigen::VectorXd&, Eigen::VectorXd&, 
	    Eigen::VectorXd&, std::function<double(double)>);



// ODE solvers
//
// first, some function wrapper typedefs
//
using ode_derivs =
  std::function<void(double, Eigen::VectorXd&, Eigen::VectorXd&)>;

using ode_integrator =
  std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, int,
		     double&, double, double, Eigen::VectorXd&,
		     double&, double&, ode_derivs)>;

using ode_rk =
  std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, int, double, double, 
		     Eigen::VectorXd&, ode_derivs)>;


using symp_derivs =
  std::function<void(double, Eigen::VectorXd&, Eigen::VectorXd&,
		     Eigen::VectorXd&)>;

using symp_integrator =
  std::function<void(Eigen::VectorXd&, Eigen::VectorXd&,
		     Eigen::VectorXd&, Eigen::VectorXd&, 
		     double, double, symp_derivs)>;

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
