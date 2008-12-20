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

#ifndef STANDALONE
#include <global.H>
#include "expand.h"
#include <global.H>
#else
#include <Particle.h>		// Stripped down particle structure
extern int this_step;
extern int Mstep;
extern int mstep;
extern unsigned multistep;
extern vector<int> stepL, stepN;
extern pthread_mutex_t coef_lock;
#endif

//! Encapsulatates a SLGridSph (Sturm-Liouville basis) for use as force method
class EmpCylSL
{
private:
  int NMAX;
  int LMAX;
  int MMAX;
  int NORDER;
  int NKEEP;

  int hallfreq, hallcount;
  string hallfile;

  double YMIN, YMAX;
  double dX, dY;
  int M, cylused, cylused1;
  vector<double> cylmass1,  cylmassE;
  bool cylmass_made;
  double cylmass;

  vector<double> r, d, m, p;

  double ASCALE;
  double HSCALE;
  double pfac, dfac, ffac;

  Matrix *facC, *facS;

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

  vector<Vector**> accum_cos0, accum_sin0;
  vector<Vector**> accum_cosL, accum_cosN;
  vector<Vector**> accum_sinL, accum_sinN;
  vector< vector<unsigned> > howmany1;
  vector<unsigned> howmany;
  vector<unsigned> dstepL, dstepN;

  Vector* accum_cos;
  Vector* accum_sin;

  Vector** accum_cos2;
  Vector** accum_sin2;
  Matrix *vc, *vs;

  Matrix tabp, tabf, tabd;

  Vector* hold;

  vector<short> coefs_made;
  bool eof_made;

  SphericalModelTable* make_sl();
  SphericalModelTable* model;

  void make_grid();
  void send_eof_grid();
  void receive_eof(int request_id, int m);
  void compute_eof_grid(int request_id, int m);
  void setup_eof_grid(void);
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

  //! No extrapolating beyond grid (default: false)
  static bool enforce_limits;

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

  //! Selector output freq
  static int HALLFREQ;

  //! Minimum radial value for basis
  static double RMIN;

  //! Maximum radial value for basis
  static double RMAX;

  //! Name of cache file
  static string CACHEFILE;

  //! Fraction of table range for basis images (for debug)
  static double HFAC;


  //! Constructor (reset must called later)
  EmpCylSL(void);

  //! Constructor with parameters
  EmpCylSL(int numr, int lmax, int mmax, int nord,
	   double ascale, double hscale);

  //! Destructor
  ~EmpCylSL(void);

  //! Reconstruct basis with new parameters
  void reset(int numr, int lmax, int mmax, int nord,
	     double ascale, double hscale);

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
  double z_to_y(double z) { return z/(fabs(z)+DBL_MIN)*asinh(fabs(z/HSCALE)); }

  //! Compute Z from non-dimensional vertical coordinate
  double y_to_z(double y) { return HSCALE*sinh(y); }

  //! Get basis function value
  void get_all(int m, int n, double r, double z, double phi,
	       double& p, double& d, double& fr, double& fz, double& fp);

  //! Setup for accumulated coefficients
  //@{
  //! All levels
  void setup_accumulation(void);
  //! Single level
  void setup_accumulation(int mlevel);
  //! For EOF
  void setup_eof(void);
  //! Clear mass counter
  void reset_mass(void);
  //@}

  //! Make coefficients from accumulated data
  //@{
  //! All levels
  void make_coefficients(void);
  //! Single level
  void make_coefficients(int mlevel);
  //! Make empirical orthgonal functions
  void make_eof(void);
  //! True if coefficients are made at all levels
  bool coefs_made_all() 
  {
    for (unsigned M=0; M<=multistep; M++) 
      if (!coefs_made[M]) return false;
    return true;
  }
  //@}

  //! Necessary member function currently unused (change design?)
  void determine_coefficients() {};
  //! Necessary member function currently unused (change design?)
  void determine_acceleration_and_potential() {};

  //! Accumulate coefficients from particle distribution
  void accumulate(vector<Particle>& p, int mlev=0, bool verbose=false);

  //! Make EOF from particle distribution
  void accumulate_eof(vector<Particle>& p, bool verbose=false);

  //! Add single particle to coefficients
  void accumulate(double r, double z, double phi, double mass, int id, int mlev=0);

  //! Add single particle to EOF coefficients
  void accumulate_eof(double r, double z, double phi, double mass, int id, int mlev=0);


  //! Evaluate potential and force field 
  void accumulated_eval(double r, double z, double phi, double& p0,
			double& p, double& fr, double& fz, double& fp);

  //! Evaluate density field
  double accumulated_dens_eval(double r, double z, double phi, double& d0);

  /** Extrapolate and sum coefficents per multistep level to get
      a complete set of coefficients for force evaluation at an
      intermediate time step
  */
  void compute_multistep_coefficients(unsigned mlevel);

  //! For updating levels
  //@{
  vector< vector<Matrix> > differS1, differC1;
  vector<double> workC1, workC, workS1, workS;
  //@}

  /** Update the multi time step coefficient table when moving particle 
      <code>i</code> from level <code>cur</code> to level 
      <code>next</code>
  */
  //@{
  virtual void multistep_update_begin();
  virtual void multistep_update(int from, int to, double r, double z, double phi, double mass, int id);
  virtual void multistep_update_finish();
  virtual void multistep_reset();
  //@}

  //! Print debug info
  void multistep_debug();

  //! Dump out coefficients to stream
  void dump_coefs(ostream& out);

  //! Dump out coefficients to stream in bianry format
  void dump_coefs_binary_last(ostream& out, double time);

  //! Dump out coefficients to stream in bianry format
  void dump_coefs_binary_curr(ostream& out, double time);

  //! Plot basis
  void dump_basis(const string& name, int step);

  //! Plot full fields for debugging
  void dump_images(const string& OUTFILE,
		   double XYOUT, double ZOUT, int OUTR, int OUTZ,
		   bool logscale);

  //! Plot basis images for debugging
  void dump_images_basis(const string& OUTFILE,
			 double XYOUT, double ZOUT, 
			 int OUTR, int OUTZ, bool logscale,
			 int M1, int M2, int N1, int N2);

  /** @name Utility functions */
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

  //! Return current value of disk scale length
  double get_ascale(void) { return ASCALE; }

  //! Return current value of disk scale height
  double get_hscale(void) { return HSCALE; }

  //! Set frequency and file name for selector output
  inline void setHall(string file, int n=50) {
    hallfile = file;
    hallfreq = n;
  }

  vector<double> sanity() { 
    vector<double> ret;
    for (int m=0; m<=MMAX; m++) ret.push_back(accum_cos[0][m]);
    return ret;
  }

};

extern void legendre_R(int lmax, double x, Matrix& p);
extern void dlegendre_R(int lmax, double x, Matrix &p, Matrix &dp);
extern void sinecosine_R(int mmax, double phi, Vector& c, Vector& s);


#endif
