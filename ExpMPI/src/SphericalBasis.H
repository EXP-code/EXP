#ifndef _SphericalBasis_H
#define _SphericalBasis_H

#include <AxisymmetricBasis.H>

//! Defines a spherical basis with members specific to spherical symmetry
class SphericalBasis : public AxisymmetricBasis
{

protected:

  //! Will cache values of the kernel required for recursion relations
  Matrix krnl;
  //! Work vectors to make recursion evaluations thread safe
  Vector *u, *du;
  //! Matrices per thread for obtaining field, legendre and coefficients
  Matrix *potd, *dpot, *legs, *dlegs, *expcoef0;
  //! Work vectors for sines and cosines for all values <code>m</code>
  Vector *cosm, *sinm;
  //! Work matrices
  Matrix dend, factorial;

  //! Initialize method
  virtual void initialize();
  //! Provided by derived classes to check consistency
  virtual void check_range();

  //! Should be called by derived classes after initialization
  void setup();

  //! Get potential and accel field values from coefficients
  void get_pot_coefs(int l, Vector& coef, double *p, double *dp);

  /** Get potential and field values from coefficients 
      (provide workspace for thread safety) */
  void get_pot_coefs_safe(int l, Vector& coef, double *p, double *dp,
			  Matrix& potd1, Matrix& dpot1);

  //! Get density from coefficient vector
  void get_dens_coefs(int l, Vector& coef, double* p);

  //! Flag whether or not work space has been initialized
  bool firstime_coef, firstime_accel;

  //! Thread method for coefficient accumulation
  virtual void * determine_coefficients_thread(void * arg);

  //! Thread method for accerlation compuation
  virtual void * determine_acceleration_and_potential_thread(void * arg);

				// Parameters
  //! The maximum radius and any scale factor
  double rmax, scale;

  //! Flag PCA selection, self_consitency and drop l=1 term
  bool selector, self_consistent, NO_L1, EVEN_L;

public:

  //! Constructor (<param>line</param> contains the paramter list)
  SphericalBasis(string& line);

  //! Destructor
  SphericalBasis::~SphericalBasis();
  
  //! Required member to compute accelration and potential
  /** This implemenation should be fine for most cases */
  virtual void get_acceleration_and_potential(vector<Particle>* Particles);

  //! Required member to compute coeifficients with threading
  /** The thread member must be supplied by the derived class */
  virtual void determine_coefficients(void);

  //! Required member to compute accleration and potential with threading
  /** The thread member must be supplied by the derived class */
  virtual void determine_acceleration_and_potential(void);

  //! Get derivative of potential
  /** 
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param dp will be returned array in harmonics l and radial order n
    for the derivative of the potential
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_dpotl(int lmax, int nmax, double r, Matrix& p, Matrix& dp, 
		 int tid) = 0;

  //! Get derivative of potential
  /**
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_potl(int lmax, int nmax, double r, Matrix& p,
		int tid) = 0;

  //! Get derivative of potential
  /** 
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param d will be returned array in harmonics l and radial order n
    for the density
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_dens(int lmax, int nmax, double r, Matrix& d,
		int tid) = 0;

  //! Get derivative of potential
  /**
    \param lmax is the maximum harmonic order
    \param nmax is the maximum radial order
    \param r is the evaluation radius
    \param p will be returned array in harmonics l and radial order n
    for the potential
    \param d will be returned array in harmonics l and radial order n
    for the density
    \param tid is the thread enumerator that allows the function to
    provide some thread safety
  */
  virtual 
  void get_potl_dens(int lmax, int nmax, double r, Matrix& p, Matrix& d,
		     int tid) = 0;

  //! Returns the potential, its derivatives & density in spherical coords
  void determine_fields_at_point_sph(double r, double theta, double phi,
				     double *tdens, double *tpotl, 
				     double *tpotr, double *tpott, double *tpotp);


  //! Returns the potential, its derivatives & density in cylindrical coords
  void determine_fields_at_point_cyl(double r, double z, double phi,
				     double *tdens, double *tpotl, 
				     double *tpotr, double *tpotz, double *tpotp);


  //! Normalization for recursion relation basis
  virtual double norm(int, int) { return 1.0; }

  //! Kernel value for recursion relation basis
  virtual double knl(int, int) { return 1.0; }

  //! Dump current coefficients
  void dump_coefs(ostream &);

};

#endif
