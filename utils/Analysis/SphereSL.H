#ifndef SphereSL_H
#define SphereSL_H 1

#include <Vector.h>
#include <string>

#include <SLGridMP2.h>
#include <localmpi.h>

/**
   Uses SLGridSph basis to evaluate expansion coeffients and provide
   potential and density basis fields
*/
class SphereSL
{

private:

  SLGridSph *sl;
  int lmax, nmax;

  Matrix potd, dpot, dpt2, dend;
  Matrix legs, dlegs, d2legs;

  Matrix cc;			// Will hold covariance matrix
  Vector cw;			// Work vector
  int nmat;

  bool coefs_defined;
  Matrix factorial;
  Matrix expcoef;
  double rscl;
  int used;

  bool compute_covar;

  Matrix covar;
  Matrix var;

  Vector work1, work2;

  void bomb(char *s);

public:

  //! Number of radial points in table (default: 800)
  static int NUMR;

  //! Use MPI?  (default: false)
  static bool mpi;

  /** Constructor
      \param mod instance of SphericalModelTable to condition eigenfunctions
      \param LMAX is the maximum harmonic order
      \param NMAX is the maximum radial order
      \param COVAR true will evaluate the covariance matrix for coeffficients
  */
  SphereSL(SphericalModelTable* mod, int LMAX, int NMAX, bool COVAR=false);
  
  //! Destructor
  ~SphereSL(void);
  
  //! Prescaling factor
  void set_scale(const double scl) { rscl = scl; }

  //! Zero out coefficients to prepare for a new expansion
  void reset_coefs(void);

  //! Make coefficients after accumulation
  void make_coefs(void);

  //! Install coefficients
  void install_coefs(Matrix& newcoef);

  //! Retrieve coefficients
  Matrix retrieve_coefs(void) { return expcoef; }

  //! Accumulate new coefficients
  void accumulate(double x, double y, double z, double mass);
  
  //! Evaluate potential and density
  void dens_pot_eval(double r, double costh, double phi,
		     double& dens0, double& dens, 
		     double& potl0, double& potl,
		     int L1=0, int L2=10000);

  //! Evaluate potential and force
  void pot_force_eval(double r, double costh, double phi,
		      double& potl, double& potr, double& pott, double& potp,
		      int L1=0, int L2=10000);
  
  void all_eval(double r, double costh, double phi,
		double& den0, double& den1,
		double& pot0, double& pot1,
		double& potr, double& pott, double& potp,
		int L1=0, int L2=10000);

  //! Return current maximum harmonic order in expansion
  int getLmax() { return lmax; }

  //! Return current maximum order in radial expansion
  int getNmax() { return nmax; }

  // @{ Evaluate Legendre polynomials and derivatives
  void legendre_R(int lmax, double x, Matrix &p);
  void legendre_R(int lmax, double x, Matrix &p, Matrix &dp);
  void legendre_R(int lmax, double x, Matrix &p, Matrix &dp, Matrix& d2p);
  // @}

  //! Access to basis
  SLGridSph *basis(void) { return sl; }

};



#endif
