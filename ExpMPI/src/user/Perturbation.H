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

  AxiSymModel *model;
  AxiSymBiorth *biorth;
  bool computed;
  bool omega_computed;

  string ID;

  int l, m, nmax;

  int Epts, Kpts, Recs, Rnum;
  double OMPI, omega;

  Vector bcoef;
  CMatrix working, current, asymp;
  void make_response(void);

  const static int default_Epts=5;
  const static int default_Kpts=4;
  const static int default_RECS=400;
  const static double default_A=1.0e-2;
  const static double default_OMPI=0.03;

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

  //! Compute the expansion coefficients
  virtual void compute_coefficients() = 0;

  //! Compute the frequency from the background model
  virtual void compute_omega() = 0;
};

#endif /* _Perturbation_h */
