// This may look like C code, but it is really -*- C++ -*-

#ifndef _Perturbation_h

#define _Perturbation_h 1

#include <massmodel.h>
#include <biorth.h>

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

  static int NINT;
  static int L1MAX;
  static bool selfgrav;
  static bool verbose;

				// Constructor
  Perturbation(int Nmax);
  virtual ~Perturbation();

				// Setup

  void set_min_ompi(double ompi) { OMPI = ompi; }
  void set_respmat_parameters(int epts, int kpts, int recs);

				// Access

  int L(void) { return l; }
  int Nmax(void) { return nmax; }
  double Omega(void) { if (!omega_computed) compute_omega(); return omega; }


  void compute_perturbation(AxiSymModel *halo_model, 
			    AxiSymBiorth *halo_ortho,
			    CVector& total, CVector& pp);

  virtual void compute_coefficients() = 0;
  virtual void compute_omega() = 0;
};

#endif /* _Perturbation_h */
