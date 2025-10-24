#ifndef _QPDistF_H
#define _QPDistF_H

#include <iomanip>
#include <string>
#include <map>

// #include "biorth.H"
#include "massmodel.H"
#include "DiskWithHalo.H"
#include "orbit.H"

/**
   Compute distribution function for axisymmetric models 
   using quadratic programming techniques
 */
class QPDistF
{
private:

				// Model
  bool nt;
  std::shared_ptr<AxiSymModel> t;
  std::shared_ptr<SphericalOrbit> orb;

				// Parameters
  double RMMAX, REMAX;
  int EGRID, KGRID, MGRID;
  double SIGMA;
  double LAMBDA, ALPHA, BETA, GAMA;
  double ROFF, EOFF, KOFF;
  double KMIN, KMAX;
  int NINT, NUMT;

				// Class variables
  bool df_computed;
  Eigen::VectorXd sigma_E, sigma_K;	// Kernel parameters
  Eigen::VectorXd Egrid, Kgrid;		// Kernel eval. points
  Eigen::VectorXd X;			// Solution vector
  double obj0;			// Objective function sans anisotropy
				//      constraint
  double obj;			// Total objective function
  int IFAIL;			// return flag from QLD

  int verbose;			// Print out solution if set;

  int NJMAX;			// Eigen::VectorXd grid for JMax
  double Emin, Emax, TOLE;
  Eigen::VectorXd JMAXE, JMAX, JMAX2;
				// Cumulative grid
  vector<double> pdf, cdf;
  int NUME, NUMK;
  double TOLK;
  typedef pair<int, int> ip;
  typedef pair<double, ip> elem;
  multimap<double, ip> realz;

				// Class utility functions
  double kernel(double x, double y, double x0, double y0, 
		double sx, double sy);
  double kernel_x(double x, double y, double x0, double y0, 
		double sx, double sy);
  double kernel_y(double x, double y, double x0, double y0, 
		double sx, double sy);
  double kernel_x2(double x, double y, double x0, double y0, 
		double sx, double sy);
  double kernel_y2(double x, double y, double x0, double y0, 
		double sx, double sy);
  double kernel_xy(double x, double y, double x0, double y0, 
		double sx, double sy);

				// Error function

  void bomb(char *s) {
    cerr << "QPDistF: " << s << endl;
    exit(-1);
  }


public:

  static bool MassEGrid;	// Choose basis by Phi(R(M))   (Default: true)
  static bool MassLinear;	// Choose linear range in mass (Default: true)
  static int ITERMAX;		// Default: 1000
  static double ITERTOL;	// Default: 1.0e-6
  static double FSIGE;		// Energy kernal variance prefactor
				// Default: 1.0
  static double FSIGK;		// J/Jmax kernal variance prefactor
				// Default: 2.0


  /** Constructor
      @param T is a pointer to an axisymmetric model (2 or 3 dim)
      @param H is a pointer to an axisymmetric model (3-d halo)
      @param rmmax is the maximum radius for the radial QP grid
      @param remax determines the maximum energy Phi(remax) for the energy grid
      @param egrid is the number of energy knots for Gaussian kernels
      @param kgrid is the number of kappa knots for Gaussian kernels
      @param mgrid is the number of radial knots for density evaluation
      @param sigma is the fraction over the spacing for the kernel width
      @param lambda is the penalty coefficient in powers of kappa
      @param alpha is the exponent for powers of kappa
      @param beta is the power law exponent for radial-grid spacing in mass fraction
      @param gama is the power law exponent for energy-grid spacing
      @param roff is the offset of the mass zero point
      @param eoff is the offset of the energy grid from Emin
      @param koff is the offset of the energy grid from Kmin
      @param kmin is the minimum value of kappa
      @param kmax is the maximum value of kappa
      @param nint is the number of knots in the velocity magnitude integral
      @param numt is the number of knots in the velocity angular integral
   */
  QPDistF(AxiSymModPtr T, AxiSymModPtr H, double rmmax, double remax, 
	  int egrid, int kgrid, int mgrid,
	  double sigma=2.0, double lambda=0.0, double alpha=-4.0, 
	  double beta=1.0, double gama=1.0,
	  double roff=0.01, double eoff=0.5, double koff=0.5, 
	  double kmin=0.0, double kmax=1.0,
	  int nint=80, int numt=40);

  QPDistF(AxiSymModPtr T, double rmmax, double remax, 
	  int egrid, int kgrid, int mgrid,
	  double sigma=2.0, double lambda=0.0, double alpha=-4.0, 
	  double beta=1.0, double gama=1.0,
	  double roff=0.01, double eoff=0.5, double koff=0.5, 
	  double kmin=0.0, double kmax=1.0,
	  int nint=80, int numt=40);

  QPDistF(AxiSymModPtr T, string file);

  ~QPDistF(void);

  void compute_distribution(void);

  void set_verbose(void);
  void set_verbose(int verbosity);

  double distf (double E, double L);
  double dfdE  (double E, double L);
  double dfdL  (double E, double L);
  double d2fdE2(double E, double L);

  double distf_EK (double E,double K);
  double dfdE_EK  (double E,double K);
  double dfdK_EK  (double E,double K);
  double d2fdE2_EK(double E,double K);
  double d2fdK2_EK(double E,double K);
  double d2fdEK_EK(double E,double K);

  void get_objective(double* OBJ0, double* OBJ, int* IFLG)
    {*OBJ0=obj0; *OBJ=obj; *IFLG=IFAIL;}

  // Cumulative grid
  void make_cdf(int ENUM, int KNUM, double KTOL=1.0e-3);
  void dump_cdf(const string& file);
  pair<double, double> gen_EK(double r1, double r2);
  
  // Read in already computed distribution function
  void read_state(string& name);

  // Write out distribution function for future use
  void write_state(string& name);

};

#endif				// QPDistF.H
