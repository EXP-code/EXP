#include <iostream>
#include <string>

#include <orbit.h>
#include <massmodel.h>
#include <biorth.h>
#include <CVector.h>
#include <interp.h>
#include <numerical.h>

#ifndef _ResPot_H
#define _ResPot_H

/**
   Compute Resonance Potential for given expansion term and perturbation
   using a mapping
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
    vector<double> dWE;
    vector<double> dWK;
    double O1, O2, Jm, dJm, E, K, I1;
    int num;
  };
  
  void check_rw(RW&);

  typedef vector<RW> ovector;
  typedef vector<ovector> o2vector;

  o2vector orbmat;

  double Emax, Emin, Kmin, Kmax;
  bool grid_computed, second_order;
  void compute_grid();

  typedef vector<double> dvector;
  bool actions_computed;

  int numx;
  double dX;
  vector<dvector> I1X, EX;
  vector<double> EminX;


  bool getValues(double I1, double I2,
		 double& O1, double& O2);
  
  bool getValues(double I1, double I2, CVector& bcoef,
		 double& O1, double& o2,
		 double& Jm, double& dJm,
		 KComplex& Ul, KComplex& dUldE, KComplex& dUldK);
  
  double xJ(double J, double Jmin, double Jmax);
  double Jx(double x, double Jmin, double Jmax);
  double dxJ(double J, double Jmin, double Jmax);

  double Rmin, Rmax, minJ, maxJ;

  void getInterp(double I1, double I2, int& indxX, int indxE[2], double cX[2], double cE[2][2], bool& noboundary);

  int ngrid;

public:

  //! Power for ang. mom. scaling (default: 0.25)
  static double ALPHA;
  //! Fractional energy offset for grid (default: 0.01)
  static double DELTA;
  //! Fractional steps for differenced derivatives
  static double DELE, DELK;
  static double TOLITR;
  static int NUME;
  //! Number of scaled ang mom points
  static int NUMX;
  //! Number of knots for orbit table
  static int RECS;
  //! Maximum number of iterations in mapping
  static int ITMAX;
  static KComplex I;


  //! Constructor
  ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	 int l, int m, int l1, int l2, int n);

  //! Destructor
  ~ResPot();

  //! Coordinate transform: Cartesian to Action-Angle
  bool coord(double* pos, double* vel,
	     double& E, double& K, double& I1, double& J,
	     double& O1, double& O2,
	     double& W1, double& W2, double& W3, 
	     double& F, double& BETA, double& PSI);
  
  //! Coordinate transform: Action-Angle to Cartesian
  bool coord(double* pos, double* vel,
	     double I1, double I2, double beta, 
	     double w1, double w2, double w3);

  //! Update actions
  bool Update(double dt, double Phase, double Omega, 
	      double amp, CVector& bcoef,
	      double* posI, double* velI,
	      double* posO, double* velO);
  
  //! Force
  bool Force(double dt, double Phase, double Omega, 
	     double amp, CVector& bcoef,
	     double* pos, double* vel, double* acc);

  //! Set first order implicit symplectic integration
  void set_1st_order() { second_order = false; }

  //! Set second order implicit symplectic integration
  void set_2nd_order() { second_order = true; }

};

double rot_matrix(int l, int m, int n, double beta);
double Ylm01(int ll, int mm);
KComplex VeeBeta(int l, int l2, int m, double beta);
extern int respot_mpi_id();

#endif

