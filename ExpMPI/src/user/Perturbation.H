// This may look like C code, but it is really -*- C++ -*-

#ifndef _Perturbation_h

#define _Perturbation_h 1

#include <massmodel.h>
#include <biorth.h>

/**
   Abstract perturbation applied through a biorthogonal expansion
*/
class Perturbation  {

protected:

  //! Cached model instance
  AxiSymModel *model;
  //! Cached biorthogonal function instance
  AxiSymBiorth *biorth;
  //! True if perturbation has been computed
  bool computed;
  //! True if perturbation frequency has been computed
  bool omega_computed;

  //! String identifying type of perturbation
  string ID;

  //! Harmonic orders and radial terms
  int l, m, nmax;

  //! Points for numerical grids
  int Epts, Kpts, Recs, Rnum;
  //! Frequencies
  double OMPI, omega;

  //! Perturabtion vector
  Vector bcoef;
  //! Temporary storage
  CMatrix working, current, asymp;
  //! Method to compute response matrix
  void make_response(void);

  //! Default value for Energy grid
  static const int default_Epts=5;
  //! Default value for Kappa grid
  static const int default_Kpts=4;
  //! Default value for number of points in orbit table
  static const int default_RECS=400;
  //! Constant imaginary part of forcing frequency (can't make this const because of portability standard)
  static double default_OMPI;

public:

  //! Number of integration knots
  static int NINT;
  //! Maximum in "l_1" expansion
  static int L1MAX;
  //! Use self gravity 
  static bool selfgrav;
  //! Verbose comments
  static bool verbose;

				// Constructor
  //! Constructor
  /**
     @param Nmax is the radial expansion order
  */
  Perturbation(int Nmax);
  //! Destructor
  virtual ~Perturbation();

				// Setup

  //! Set constant imaginary part of frequency (for avoiding poles on real axis)
  void set_min_ompi(double ompi) { OMPI = ompi; }

  //! Set parameters for self-gravitiationg solution
  void set_respmat_parameters(int epts, int kpts, int recs);

				// Access

  //! Return harmonic order
  int L(void) { return l; }

  //! Return radial order
  int Nmax(void) { return nmax; }

  //! Return current frequency
  double Omega(void) { if (!omega_computed) compute_omega(); return omega; }

  //! Compute (recompute) the perturbation
  void compute_perturbation(AxiSymModel *halo_model, 
			    AxiSymBiorth *halo_ortho,
			    CVector& total, CVector& pp);


  //! Register halo model with perturbation
  void set_model(AxiSymModel *halo_model) { model = halo_model; }


  //! Compute the expansion coefficients
  virtual void compute_coefficients() = 0;

  //! Compute the frequency from the background model
  virtual void compute_omega() = 0;

  //! Get value of perturbation potential
  virtual double eval(double r) = 0;
};

#endif /* _Perturbation_h */
