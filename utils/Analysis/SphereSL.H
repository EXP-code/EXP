#ifndef SphereSL_H
#define SphereSL_H 1

#include <string>

#include <SLGridMP2.H>
#include <localmpi.H>
#include <Eigen/Eigen>

/**
   Uses SLGridSph basis to evaluate expansion coeffients and provide
   potential and density basis fields
*/
class SphereSL
{

private:

  std::shared_ptr<SLGridSph> sl;
  int lmax, nmax;
  double rmin, rmax;

  Eigen::MatrixXd potd, dpot, dpt2, dend;
  Eigen::MatrixXd legs, dlegs, d2legs;

  bool coefs_defined;
  Eigen::MatrixXd factorial;
  Eigen::MatrixXd expcoef;
  double rscl;
  int used;

  bool compute_covar;

  using matT = std::vector<Eigen::MatrixXd>;
  using vecT = std::vector<Eigen::VectorXd>;

  matT covar;
  vecT mean;
  vecT svar;
  matT uvec;
  double totalMass;

  std::vector<vecT> meanB;
  std::vector<double> massB;
  int npart, curbin;

  double minSNR, maxSNR;

  Eigen::VectorXd work;

  void bomb(char *s);

  //! Coefficient magic number
  const unsigned int cmagic = 0xc0a57a2;

  //! For coefficient writing
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
  EigenColMajor;

public:

  //! Number of eigenfunctions for each subspace for PNG figure
  static int NEV;

  //! Number of radial points in table (default: 800)
  static int NUMR;

  //! Use MPI?  (default: false)
  static bool mpi;

  //! Hall exponent (default: 1)
  static double HEXP;

  /** Constructor
      \param mod instance of SphericalModelTable to condition eigenfunctions
      \param LMAX is the maximum harmonic order
      \param NMAX is the maximum radial order
      \param CMAP is the SLGridSph coordindate mapping type
      \param SCALE is the SLGridSph coordindate scale length
      \param COVAR true will evaluate the covariance matrix for coeffficients
      \param NPART >0 will evaluate the covariance matrix using a jackknife partition with number given
  */
  SphereSL(std::shared_ptr<SphericalModelTable> mod, int LMAX, int NMAX,
	   int CMAP=0, double SCALE=0.067, bool COVAR=false, int NPART=0);
  
  //! Destructor
  ~SphereSL(void);
  
  //! Prescaling factor
  void set_scale(const double scl) { rscl = scl; }

  //! Zero out coefficients to prepare for a new expansion
  void reset_coefs(void);

  //! Make coefficients after accumulation
  void make_coefs(void);

  //! Covariance analysis
  void make_covar(bool verbose=false);

  //! Install coefficients
  void install_coefs(Eigen::MatrixXd& newcoef);

  //! Retrieve coefficients
  Eigen::MatrixXd retrieve_coefs(void) { return expcoef; }

  //! Minimum SNR coefficient value
  double getMinSNR(void) { return minSNR; }

  //! Maximum SNR coefficient value
  double getMaxSNR(void) { return maxSNR; }

  //! Get mass on grid
  double getMass(void)   { return totalMass; }

  //! Retrieve coefficients with SNR trim
  Eigen::MatrixXd get_trimmed(double snr, double mass, bool Hall=false);

  //! Get fractional power in trimmed coefficients
  double get_power(double snr, double mass);

  //! Retrieve coefficients
  std::vector<Eigen::MatrixXd> retrieve_covar(void) { return covar; }

  //! Retrieve potential 
  Eigen::MatrixXd get_pot(double x) { sl->get_pot(potd,x); return potd;}
  //sl->get_dens (dend, r/rscl);
  //sl->get_pot  (potd, r/rscl);
  //sl->get_force(dpot, r/rscl);


  //! Accumulate new coefficients
  void accumulate(double x, double y, double z, double mass);
  
  //! Evaluate potential and density
  void dens_pot_eval(double r, double costh, double phi,
		     double& dens0, double& dens, 
		     double& potl0, double& potl,
		     int L1=0, int L2=10000, int N1=0, int N2=10000);

  //! Evaluate potential and force
  void pot_force_eval(double r, double costh, double phi,
		      double& potl, double& potr, double& pott, double& potp,
		      int L1=0, int L2=10000, int N1=0, int N2=10000);
  
  void all_eval(double r, double costh, double phi,
		double& den0, double& den1,
		double& pot0, double& pot1,
		double& potr, double& pott, double& potp,
		int L1=0, int L2=10000, int N1=0, int N2=10000);

  //! Return current maximum harmonic order in expansion
  int getLmax() { return lmax; }

  //! Return current maximum order in radial expansion
  int getNmax() { return nmax; }

  // @{ Evaluate Legendre polynomials and derivatives
  void legendre_R(int lmax, double x, Eigen::MatrixXd &p);
  void legendre_R(int lmax, double x, Eigen::MatrixXd &p, Eigen::MatrixXd &dp);
  void legendre_R(int lmax, double x, Eigen::MatrixXd &p, Eigen::MatrixXd &dp,
		  Eigen::MatrixXd& d2p);
  // @}

  //! Access to basis
  std::shared_ptr<SLGridSph> basis(void) { return sl; }

  //! Dump coefficients using the standard header
  void dump_coefs(double time, std::ostream& out);

};



#endif
