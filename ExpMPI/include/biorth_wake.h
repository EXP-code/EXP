// This may look like C code, but it is really -*- C++ -*-

// Wake generation class

#ifndef _biorth_wake_h
#define _biorth_wake_h 1

static const char rcsid_biorth_wake[] = "$Id$";

#include <biorth.h>
#include <CVector.h>

class BiorthWake
{

private:

  AxiSymBiorth *bio;
  int lmax, nmax;

  Matrix norm_grid;		// Stored as the sqr root of normalization
				// for efficiency in wake reconstruction

  bool coefs_defined;
  Matrix factorial;
  Matrix expcoef;
  double rscl;
  int used;

  void accumulate_2d(double x, double y, double mass);
  void accumulate_3d(double x, double y, double z, double mass);

  void reconstruct_2d(double r, double phi,
		      double& dens0, double& dens, 
		      double& potl0, double& potl,
		      int L1=0, int L2=10000);

  void reconstruct_3d(double r, double costh, double phi,
		      double& dens0, double& dens, 
		      double& potl0, double& potl,
		      int L1=0, int L2=10000);

  void bomb(const char *s);

				// Orientation

  Complex I;
  CVector ylm;
  bool init_orientation;
  double *param;
  double **ambp, *amby, *ptry, *psum;

  int iter, ll, mm;
  double tol;
  static int iterDef;
  static double tolDef;
  static int ndim;

  int nfunk;
  void modulo_param(double *params);
  void orientation_init(void);
  double amoeba_energy(double *params);
  void get_transform(double& phi, double& theta, double& psi, double& cost);
  void amoeba(void);
  double amotry(int ihi, double fac);
				// Debugging routines only
  void test_transform(void);
  void check_orientation(double, double, double);
  Complex test_fct(double, double);

public:

// Constructor
  BiorthWake(AxiSymBiorth *BIO, int LMAX, int NMAX);

// Destructor
  ~BiorthWake(void);

// Member functions
  void set_scale(const double scl) { rscl = scl; }
  void reset_coefs(void);

  void accumulate(double x, double y, double z, double mass);

  void reconstruct(double r, double costh, double phi,
			     double& dens0, double& dens, 
			     double& potl0, double& potl,
			     int L1=0, int L2=10000);

  void orientation(int L, int M, Vector& phi, Vector& theta, Vector& psi,
		   Vector& cost);

  int get_amoeba_iterations(void) { return nfunk; }

  double energy(double *params);
};


#endif
