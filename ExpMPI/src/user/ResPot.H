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

  LegeQuad* lq;

public:

  static bool use_grid;
  static double DELTA;
  static int NREC;
  static int NUME;
  static int NUMK;
  static int NUMT;
  static KComplex I;


  //! Constructor
  ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	 int l, int m, int l1, int l2, int n);

  //! Destructor
  ~ResPot();

  //! Get potential
  double Pot(double* pos, double* vel, double Omega, double time,
	     CVector &bcoef);

};

double rot_matrix(int l, int m, int n, double beta);
double Ylm01(int ll, int mm);
KComplex VeeBeta(int l, int l2, int m, double beta);

#endif

