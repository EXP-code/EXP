#include <iostream>
#include <string>

#include <orbit.h>
#include <massmodel.h>
#include <biorth.h>
#include <CVector.h>
#include <interp.h>

#ifndef _ResPot_H
#define _ResPot_H

class EKinterp : public Matrix
{
private:
  
  double emin, emax, kmin, kmax, dE, dK;
  int nume, numk;
  
public:

  //! Constructor
  EKinterp();
  EKinterp(double Emin, double Emax, double Kmin, double Kmax,
	   int numE, int numK);

  //! Set or reset array
  void reset(double Emin, double Emax, double Kmin, double Kmax,
	   int numE, int numK);

  //! Interpolate
  double value(double E, double K);
  double deriv(double E, double K);
  void eval(double E, double K, double& value, double& deriv);
};


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
  bool grid_computed;
  void compute_grid();

  typedef vector<double> dvector;
  bool actions_computed;
  double I1min, I1max;
  EKinterp I1_EK, O1_EK, O2_EK;
  dvector I1X;
  vector<dvector> I1E, I1J;

  void compute_actions();
  double getE(double I1, double I2);
  void getPert(double E, double K, CVector& bcoef,
	       double& Jm, double& dJm,
	       KComplex& Ul, KComplex& dUldE, KComplex& dUldK);

  double xJ(double J, double Jmin, double Jmax);
  double Jx(double x, double Jmin, double Jmax);
  double dxJ(double J, double Jmin, double Jmax);

  double dE, dK, Rmin, Rmax;
  vector<double> EE, KK, Jmax;

  int ngrid;

public:

  static double ALPHA;
  static double DELTA;
  static int NREC;
  static int NUME;
  static int NUMK;
  static int NUMT;
  static int NUMI;
  static int ITMAX;
  static KComplex I;


  //! Constructor
  ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	 int l, int m, int l1, int l2, int n);

  //! Destructor
  ~ResPot();

  //! Coordinate transform: Cartesian to Action-Angle
  int coord(double* pos, double* vel,
	    double& E, double& K, double& J,
	    double& W1, double& W2, double& W3, 
	    double& F, double& BETA, double& PSI);

  //! Coordinate transform: Action-Angle to Cartesian
  int coord(double* pos, double* vel,
	    double I1, double I2, double beta, 
	    double w1, double w2, double w3);

  //! Update actions
  int Update(double dt, double Phase, double Omega, 
	     double amp, CVector& bcoef,
	     double* posI, double* velI,
	     double* posO, double* velO);

  //! Force
  int Force(double dt, double Phase, double Omega, 
	    double amp, CVector& bcoef,
	    double* pos, double* vel, double* acc);

};

double rot_matrix(int l, int m, int n, double beta);
double Ylm01(int ll, int mm);
KComplex VeeBeta(int l, int l2, int m, double beta);

#endif

