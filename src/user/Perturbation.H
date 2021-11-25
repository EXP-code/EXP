// This may look like C code, but it is really -*- C++ -*-

#ifndef _Perturbation_h

#define _Perturbation_h 1

#include <complex>
#include <Eigen/Eigen>

#include <massmodel.H>
#include <biorth.H>

/**
   Abstract perturbation applied through a biorthogonal expansion
*/
class Perturbation  {

protected:

  //! Cached model instance
  std::shared_ptr<AxiSymModel> model;

  //! Cached biorthogonal function instance
  std::shared_ptr<AxiSymBiorth> biorth;

  //! True if perturbation has been computed
  bool computed;

  //! True if perturbation frequency has been computed
  bool omega_computed;

  //! String identifying type of perturbation
  std::string ID;

  //! Harmonic orders and radial terms
  //@{
  int l, m, nmax;
  //@}

  //! Points for numerical grids
  //@{
  int Epts, Kpts, Recs, Rnum;
  //@}

  //! Frequencies
  //@{
  double OMPI, omega;
  //@}

  //! Perturabtion vector
  Eigen::VectorXd bcoef;

  //! Temporary storage
  //@{
  Eigen::MatrixXcd working, current, asymp;
  //@}
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
  /** Constructor
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
  void compute_perturbation(AxiSymModPtr halo_model, 
			    AxiSymBioPtr halo_ortho,
			    Eigen::VectorXcd& total, Eigen::VectorXcd& pp);


  //! Register halo model with perturbation
  void set_model(AxiSymModPtr halo_model) { model = halo_model; }


  //! Compute the expansion coefficients
  virtual void compute_coefficients() = 0;

  //! Compute the frequency from the background model
  virtual void compute_omega() = 0;

  //! Get value of perturbation potential
  virtual double eval(double r) = 0;
};

#endif /* _Perturbation_h */
