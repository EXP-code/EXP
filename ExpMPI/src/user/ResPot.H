#include <iostream>
#include <string>

#include <orbit.h>
#include <massmodel.h>
#include <biorth.h>
#include <gaussQ.h>
#include <CVector.h>
#include <interp.h>

#ifndef _ResPot_H
#define _ResPot_H

/**
   Computes the resonance potential for a given resonance
*/
class ResPot
{
private:
  int L, M, L1, L2, NMAX;

  AxiSymModel *halo_model;
  AxiSymBiorth *halo_ortho;

  SphericalOrbit *orb;

  class RW {
  public:
    vector<double> r;
    vector<double> w1;
    vector<double> f;
    vector<double> W;
    int num;
  };
  
  typedef vector<RW> ovector;
  typedef vector<ovector> o2vector;

  o2vector orbmat;

  double Emax, Emin, Kmin, Kmax;
  void compute_grid();

  double dE, dK;
  vector<double> EE, KK, Jmax;

  //! Numerical constant
  static KComplex I;
  
public:

  //! Relative radial offset for numerical differential (default: 0.01)
  static double DELTA;
  //! Number of ENERGY grid points for tabled grid (default: 400)
  static int NUME;
  //! Number of KAPPA grid points for tabled grid (default: 100)
  static int NUMK;

  //! Constructor
  ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	 int l, int m, int l1, int l2, int n);

  //! Destructor
  ~ResPot();

  //! Get potential
  double Pot(double* pos, double* vel, double phase, CVector &bcoef);

};

double rot_matrix(int l, int m, int n, double beta);
double Ylm01(int ll, int mm);
KComplex VeeBeta(int l, int l2, int m, double beta);

#endif

