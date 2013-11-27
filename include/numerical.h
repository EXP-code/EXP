#ifndef _numerical_h
#define _numerical_h

/* 
	some simple macros
*/

#ifndef MAX
#define MAX(bob, joe) (((bob) > (joe)) ? (bob) : (joe))
#endif
#ifndef MIN
#define MIN(bob, joe) (((bob) < (joe)) ? (bob) : (joe))
#endif

#define FLOAT(bob)  ((float) (bob))
#define DOUBLE(bob) ((double) (bob))
#define SQR(bob) ((bob)*(bob))
#define Pi 3.141592653589793
#ifndef TINY
#define TINY 1.e-10
#endif
#define Sign(bob) (((bob)>0) ? (1.0) : (-1.0))

/* special functions */

// double BesselJ(int, double);
// double BesselI(int, double);
// double BesselK(int, double);
// double gammln(double);


double poly(double, int, double *);

/* some numerical recipes utility functions */

double *nr_vector(int, int);
int *nr_ivector(int, int);
double **nr_matrix(int, int, int, int);
void free_nr_vector(double *, int, int);
void free_nr_ivector(int *, int, int);
void free_nr_matrix(double **, int, int, int, int);
void nrerror(const char *);
void sort(int, double *);



/* typedefs for 1D functions used in quadrature, minimization, and roots */

typedef double (*func_1d)(double);
typedef double (*func_2d)(double, double);
typedef void (*dfunc_1d)(double, double *, double *);





/* numerical quadrature */

void gauleg(double, double, double *, double *, int);
double qadapt(double, double, func_1d, double);
double qadapt2d(double, double, double, double, func_2d, double);



/* root finding */


double zbrent(func_1d, double, double, double);

double rtbis(func_1d, double, double, double);

double rtsafe(dfunc_1d,  double, double, double);

int zbrac(func_1d, double *, double *);




/* 1D minimization */

double brent(double, double, double, func_1d, double, double *);
void mnbrak(double *, double *, double *, double *, double *, 
	double *, func_1d);



/* random numbers */

/*
void Init_Random (int);
double Random (double, double);
double Gaussian_Variable (double);
void Random_Vector (double, double *x);
double Random_Polar_Angle (double, double);
*/



/* ODE solvers */


/* first, some typedefs for passing function pointers */

	
typedef void (*ode_derivs)(double, double *, double *);

typedef void (*ode_integrator)(double *, double *, int, double *, double, 
	double, double *, double *, double *, ode_derivs);

typedef void (*ode_rk)(double *, double *, int, double, double, 
	double *, ode_derivs);


typedef void (*symp_derivs)(double, double *, double *, double *);

typedef void (*symp_integrator)(double *, double *, double *, double *, 
				double, double, symp_derivs);


/* prototypes for ODE solvers */

void integrate_ode(double *, double, double, double *, double, int,
		   ode_derivs, ode_integrator);

void onestep(double *, double *, int, double *, double, 
		ode_derivs, ode_integrator);
	
void rkqc(double *, double *, int , double *, double, double, 
		double *, double *, double *, ode_derivs);
		
void rk2qc (double *, double *, int , double *, double, double, 
		double *, double *, double *, ode_derivs);

void bsstep (double *, double *, int , double *, 
	double, double, double *, double *, double *, ode_derivs);

void constant_step (double *, double *, int, double, ode_derivs, ode_rk);

void rk4(double *, double *, int, double, double, double *, ode_derivs);

void sia4(double *, double *, double *, double *, double, double, symp_derivs);




/* interpolation stuff */


void locate (double *, int, double, int *);

/*
#ifdef __cplusplus
	class Matrix;
	class Vector;

	Matrix build_cubic_table(Vector &, Vector &, Vector &);
	double cubic_value(Matrix &, Vector &, double, double *);

#endif
*/



/* command line parser utilities */

void init_cmdparse ( int, char **);
char *strarg (char *, int);
int iarg (char *, int);
int cmdflag_set (char *);
double farg (char *, int);

#endif
