#include <iostream>
#include <string>

#include <orbit.H>
#include <massmodel.H>
#include <biorth.H>
#include <interp.H>
#include <numerical.H>

#ifndef _ResPotOrb_H
#define _ResPotOrb_H

/**
   Compute Resonance Potential for given expansion term and perturbation
   using a mapping
*/
class ResPotOrb
{
public:
  
  //! Number of return codes
  const static int NumDesc = 11;
  
  //! Return code name (character string)
  const static char* ReturnDesc[];

  //! Return codes
  enum ReturnCode {OK, 
		   CoordRad, CoordVel, CoordCirc, CoordTP, CoordRange,
		   CoordBadPos, CoordBadVel, 
		   UpdateBadL, UpdateIterate, UpdateBadVal};
  
  
private:
  int L, M, L1, L2;

  AxiSymModPtr halo_model;
  std::shared_ptr<SphericalOrbit> orb;
  
  class RW {
  public:
    std::vector<double> r, ff, w1, f;
    std::vector<double> r_Ep, ff_Ep;
    std::vector<double> r_Em, ff_Em;
    std::vector<double> r_Kp, ff_Kp;
    std::vector<double> r_Km, ff_Km;
    double O1, O2, Jm, dJm, E, K, I1;
    int num;
  };
  
  void check_rw(RW&);
  
  typedef std::vector<RW> ovector;
  typedef std::vector<ovector> o2vector;
  
  o2vector orbmat;
  
  double Emax, Emin, Kmin, Kmax, Kupd, delE, delK;
  bool grid_computed, second_order;
  void compute_grid();

  double W, dWE, dWK;

  typedef std::vector<double> dvector;
  bool actions_computed;
  
  int numx;
  double dX;
  std::vector<dvector> I1X, EX;
  std::vector<double> EminX;
  
  
  bool getValues(double I1, double I2,
		 double& O1, double& O2);
  
  bool getValues(double rsat, double I1, double I2,
		 double& O1, double& O2,
		 double& Jm, double& dJm,
		 std::complex<double>& Ul, std::complex<double>& dUldE, std::complex<double>& dUldK);
  
  double xJ(double J, double Jmin, double Jmax);
  double Jx(double x, double Jmin, double Jmax);
  double dxJ(double J, double Jmin, double Jmax);
  
  double Rmin, Rmax, minJ, maxJ;
  double MASS;
  
  void getInterp(double I1, double I2, int& indxX, int& indxE, double cX[2], double cE[2], bool& noboundary);

  ReturnCode Update2(double dt, std::vector<double>& Phase,
		     double rsat, double amp,
		     double* posI, double* velI,
		     double* posO, double* velO);
  
  ReturnCode Update3(double dt, std::vector<double>& Phase,
		     double rsat, double amp,
		     double* posI, double* velI,
		     double* posO, double* velO);
  
  int ngrid;
  
  string dbgFILE;

public:
  
  //! Power for ang. mom. scaling (default: 0.25)
  static double ALPHA;

  //! Fractional energy offset for grid (default: 0.001)
  static double DELTA_E;

  //! Kappa offset for grid (default: 0.001)
  static double DELTA_K;

  //! Beta offset for grid (default: 0.001)
  static double DELTA_B;

  //! Fractional steps for differenced derivatives
  static double DELE;

  //! Fractional steps for differenced derivatives
  static double DELK;

  //! Maximum number of iterations
  static double TOLITR;

  //! Number of points in energy grid
  static int NUME;

  //! Number of scaled ang mom points
  static int NUMX;

  //! Number of knots for orbit table
  static int RECS;

  //! Maximum number of iterations in mapping
  static int ITMAX;

  //! Complex sqrt(-1)
  static std::complex<double> I;
  
  //! Constructor
  ResPotOrb(AxiSymModPtr mod, double mass, int l, int m, int l1, int l2);

  //! Destructor
  ~ResPotOrb();
  
  //! Coordinate transform: Cartesian to Action-Angle
  ReturnCode coord(double* pos, double* vel,
		   double& E, double& K, double& I1, double& J,
		   double& O1, double& O2,
		   double& W1, double& W2, double& W3, 
		   double& F, double& BETA, double& PSI);
  
  //! Coordinate transform: Action-Angle to Cartesian
  ReturnCode coord(double* pos, double* vel,
		   double I1, double I2, double beta, 
		   double w1, double w2, double w3);
  
  //! Update actions
  ReturnCode Update(double dt, std::vector<double>& Phase,
		    double rsat, double amp,
		    double* posI, double* velI,
		    double* posO, double* velO);
  
  //! K value from last update
  double K() { return Kupd; }
  
  //! Force
  ReturnCode Force(double dt, std::vector<double>& Phase,
		   double rsat, double amp, 
		   double* pos, double* vel, double* acc);
  
  //! Set first order implicit symplectic integration
  void set_1st_order() { second_order = false; }
  
  //! Set second order implicit symplectic integration
  void set_2nd_order() { second_order = true; }
  
  //! Set debug file
  void set_debug_file(const string& dbg) { dbgFILE = dbg; }

};

double rot_matrix(int l, int m, int n, double beta);
double Ylm01(int ll, int mm);
std::complex<double> VeeBeta(int l, int l2, int m, double beta);
extern string respotorb_mpi_id();

#endif

