// This may look like C code, but it is really -*- C++ -*-

#ifndef _EmpOrth_h
#define _EmpOrth_h

#include <vector>

#include <gaussQ.h>
#include <math.h>
#include <values.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <SLGridMP2.h>

#include "expand.h"

//! Encapsulatates a SLGridSph (Sturm-Liouville basis) for use as force method
class EmpCylSL
{
private:
  int NMAX;
  int LMAX;
  int MMAX;
  int NORDER;
  int NKEEP;

  double YMIN, YMAX;
  double dX, dY;
  int M, cylused, cylused1;
  double cylmass;

  vector<double> r, d, m, p;

  Matrix *facC, *facS;
  Vector** accum_cos0;
  Vector** accum_sin0;
  
  int rank2, rank3;

  double *MPIin, *MPIout;
  double *MPIin_eof, *MPIout_eof;
  double *mpi_double_buf1, *mpi_double_buf2, *mpi_double_buf3;
  int MPIbufsz, MPItable;
  bool MPIset, MPIset_eof;
  MPI_Status status;

  double**** SC;
  double**** SS;

  Matrix *var;

  Vector ev;
  Matrix ef;
  Matrix potd, dpot, dend;
  Vector *cosm, *sinm;
  Matrix *legs, *dlegs;

  SLGridSph *ortho;

  double Rtable, XMIN, XMAX;

  Matrix** potC;
  Matrix** densC;
  Matrix** rforceC;
  Matrix** zforceC;

  Matrix** potS;
  Matrix** densS;
  Matrix** rforceS;
  Matrix** zforceS;

  Matrix* table;

  Matrix* tpot;
  Matrix* tdens;
  Matrix* trforce;
  Matrix* tzforce;

  Vector** accum_cos;
  Vector** accum_sin;
  Vector** accum_cos2;
  Vector** accum_sin2;
  Matrix *vc, *vs;

  Matrix tabp, tabf, tabd;

  Vector* hold;

  bool coefs_made;
  bool eof_recompute;
  bool eof_made;

  SphericalModelTable* make_sl();
  SphericalModelTable* model;

  void make_grid();
  void send_eof_grid();
  void receive_eof(int request_id, int m);
  void compute_eof_grid(int request_id, int m);
  void setup_eof_grid(void);
  void setup_eof(void);
  void accumulate_eof(double r, double z, double phi, double mass, int id);
  void make_eof(void);
				// 1=write, 0=read
				// return: 0=failure
  int cache_grid(int);		
  double integral(int, int, int, int);
  void get_pot(Matrix&, Matrix&, double, double);
  void pca_hall(void);
  double massR(double R);
  double densR(double R);

  void bomb(string oops);

  pthread_mutex_t used_lock, cos_coef_lock, sin_coef_lock;

public:

  //! Type of density model to use
  enum EmpModel {
    Exponential,
    Gaussian, 
    Plummer
  };

  //! TRUE if density is computed
  static bool DENS;

  //! TRUE if signal-to-noise methods are on
  static bool SELECT;

  //! TRUE if we are using coordinate mapping
  static bool CMAP;

  //! TRUE if mapping is logarithmic
  static bool logarithmic;

  //! Density model type
  static EmpModel mtype;
  
  //! Radial basis grid in radial direction
  static int NUMX;
 
  //! Radial basis grid in vertical direction
  static int NUMY;

  //! Number of bases to print (for debug)
  static int NOUT;

  //! Number of entries in radial basis table
  static int NUMR;

  //! Minimum radial value for basis
  static double RMIN;

  //! Maximum radial value for basis
  static double RMAX;

  //! Radial scale length
  static double ASCALE;

  //! Vertical scale height
  static double HSCALE;

  //! Name of cache file
  static string CACHEFILE;

  //! Name of cache table file
  static string TABLEFILE;


  //! Constructor (reset must called later)
  EmpCylSL(void);

  //! Constructor with parameters
  EmpCylSL(int numr, int lmax, int mmax, int nord);

  //! Destructor
  ~EmpCylSL(void);

  //! Reconstruct basis with new parameters
  void reset(int numr, int lmax, int mmax, int nord);

  //! Read basis from cache file
  int read_cache(void);

  //! Parameter access: get norder
  int get_order(void) {return NORDER;}

				// Z coordinate transformation

  /*
  inline double z_to_y(double z) { return asinh(z/HSCALE); }
  inline double y_to_z(double y) { return HSCALE*sinh(y); }
  */

  //! Compute non-dimensional vertical coordinate from Z
  double z_to_y(double z) { return z/(fabs(z)+DBL_MIN)*asinh(fabs(z)/HSCALE); }

  //! Compute Z from non-dimensional vertical coordinate
  double y_to_z(double y) { return HSCALE*sinh(y); }

  //! Compute new orthogonal basis from phase space on next step
  void compute_eof(void) { eof_recompute = true; }

  //! Get basis function value
  void get_all(int m, int n, double r, double z, double phi,
	       double& p, double& d, double& fr, double& fz, double& fp);

  //! Setup for accumulated coefficients
  void setup_accumulation(void);

  //! Make coefficients from accumulated data
  void make_coefficients(void);

  //! Necessary member function currently unused (change design?)
  void determine_coefficients() {};
  //! Necessary member function currently unused (change design?)
  void determine_acceleration_and_potential() {};

  //! Accumulate coefficients from particle distribution
  void accumulate(vector<Particle>& p);

  //! Add single particle to coefficients
  void accumulate(double r, double z, double phi, double mass, int id);

  //! Evaluate potential and force field 
  void accumulated_eval(double r, double z, double phi,
			double& p, double& fr, double& fz, double& fp);

  //! Evaluate density field
  double accumulated_dens_eval(double r, double z, double phi);

  //! Dump out coefficients to stream
  void dump_coefs(ostream& out);

  //! Dump out coefficients to stream in bianry format
  void dump_coefs_binary_last(ostream& out, double time);

  //! Dump out coefficients to stream in bianry format
  void dump_coefs_binary_curr(ostream& out, double time);

  //! Plot basis
  void dump_basis(const string& name, int step);

  //! Utility
  // @{

  //! Compute Associated Legendre Polynomials, return Matrix type
  void legendre_R(int lmax, double x, Matrix& p);
  /** Compute Associated Legendre Polynomials and derivitives, 
      return Matrix type */
  void dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp);
  //! Compute vectors of sines and cosines by recursion
  void sinecosine_R(int mmax, double phi, Vector& c, Vector& s);

  // @}

  //! Convert from non-dimensional to dimensional radial coordinate
  double xi_to_r(double);

  //! Convert from dimension to non-dimension radial coordinate
  double r_to_xi(double);

  //! Jacobian
  double d_xi_to_r(double);

};

extern void legendre_R(int lmax, double x, Matrix& p);
extern void dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp);
extern void sinecosine_R(int mmax, double phi, Vector& c, Vector& s);


#endif
