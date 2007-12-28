// This may look like C code, but it is really -*- C++ -*-

//
// Compute distribution function for axisymmetric models 
// using quadratic programming techniques
//

#ifndef _QPDistF_h
#define _QPDistF_h 1

#include <biorth.h>
#include <massmodel.h>
#include <orbit.h>

class QPDistF
{
private:

				// Model
  AxiSymModel *t;
  SphericalOrbit *orb;

				// Parameters
  double RMMAX, REMAX;
  int EGRID, KGRID, MGRID;
  double SIGMA;
  double LAMBDA, ALPHA, BETA, GAMA;
  double ROFF, EOFF, KOFF;
  double KMIN, KMAX;
  int NINT, NUMT;

				// Class variables
  bool df_computed;
  Vector sigma_E, sigma_K;	// Kernel parameters
  Vector Egrid, Kgrid;		// Kernel eval. points
  Vector X;			// Solution vector
  double obj0;			// Objective function sans anisotropy
				//      constraint
  double obj;			// Total objective function
  int IFAIL;			// return flag from QLD

  int VERBOSE;			// Print out solution if set;

  int NJMAX;			// Vector grid for JMax
  double Emin, Emax, TOLE;
  Vector JMAXE, JMAX, JMAX2;
				// Class utility functions
  void compute_distribution(void);
  double kernel(double x, double y, double x0, double y0, 
		double sx, double sy);
  double kernel_x(double x, double y, double x0, double y0, 
		double sx, double sy);
  double kernel_y(double x, double y, double x0, double y0, 
		double sx, double sy);
  double* convert(Matrix& a, int mm=0, int nn=0);
  double* convert(Vector& a, int nn=0);

				// Error function

  void bomb(const char *s) {
    cerr << "QPDistF: " << s << endl;
    exit(-1);
  }


public:

  static bool MassEGrid;	// Choose basis by Phi(R(M))
				// Default: TRUE
  static int ITERMAX;		// Default: 1000
  static double ITERTOL;	// Default: 1.0e-6
  static double FSIGE;		// Energy kernal variance prefactor
				// Default: 1.0
  static double FSIGK;		// J/Jmax kernal variance prefactor
				// Default: 2.0


  QPDistF(AxiSymModel *T, double rmmax, double remax, 
	  int egrid, int kgrid, int mgrid,
	  double sigma=2.0, double lambda=0.0, double alpha=-4.0, 
	  double beta=1.0, double gama=1.0,
	  double roff=0.01, double eoff=0.5, double koff=0.5, 
	  double kmin=0.0, double kmax=1.0,
	  int nint=20, int numt=20);

  QPDistF(AxiSymModel *T, string file);

  ~QPDistF(void);

  void set_verbose(void);

  double distf(double E, double L);
  double dfdE(double E, double L);
  double dfdL(double E, double L);

  double distf_EK(double E,double K);
  double dfdE_EK(double E,double K);
  double dfdK_EK(double E,double K);

  void get_objective(double* OBJ0, double* OBJ, int* IFLG)
    {*OBJ0=obj0; *OBJ=obj; *IFLG=IFAIL;}

  // Read in already computed distribution function
  void read_state(string& name);

  // Write out distribution function for future use
  void write_state(string& name);

};

#endif				// QPDistF.h
