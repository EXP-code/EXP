#ifndef _BasisClient_H
#define _BasisClient_H

#include <SLGridMP2.H>
#include <localmpi.H>
#include <Eigen/Eigen>
#include <yaml-cpp/yaml.h>

//! Defines a basis-based potential and acceleration class
class BasicBasis
{
protected:

  //! Contains parameter database
  const YAML::Node& conf;

  /** @name Utility functions */
  // @{

  // @{ Evaluate Legendre polynomials and derivatives
  void legendre_R(int lmax, double x, Eigen::MatrixXd &p);
  void legendre_R(int lmax, double x, Eigen::MatrixXd &p, Eigen::MatrixXd &dp);
  void legendre_R(int lmax, double x, Eigen::MatrixXd &p, Eigen::MatrixXd &dp,
		  Eigen::MatrixXd& d2p);
  // @}

  //! Compute vectors of sines and cosines by recursion
  void sinecosine_R(int mmax, double phi, Eigen::VectorXd& c, Eigen::VectorXd& s);
  // @}

  //! Total mass on grid
  double totalMass;

public:

  //! Constructor
  BasicBasis(const YAML::Node& conf);

  //! Get potential, density and forces at a point
  virtual void getFields
  (double x, double y, double z,
   double& tdens0, double& tpotl0, double& tdens, double& tpotl, 
   double& tpotx, double& tpoty, double& tpotz) = 0;

  //! Accumulate new coefficients
  void accumulate(double x, double y, double z, double mass);

  //! Get mass on grid
  double getMass(void)   { return totalMass; }

  //! Zero out coefficients to prepare for a new expansion
  virtual void reset_coefs(void) = 0;

  //! Make coefficients after accumulation
  virtual void make_coefs(void) = 0;
};

/**
   Uses SLGridSph basis to evaluate expansion coeffients and provide
   potential and density basis fields
*/
class SphericalSL : public BasicBasis
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

  Eigen::VectorXd work;

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

  //! Constructor
  SphericalSL(const YAML::Node& conf);
  
  //! Destructor
  ~SphericalSL(void);
  
  //! Prescaling factor
  void set_scale(const double scl) { rscl = scl; }

  //! Zero out coefficients to prepare for a new expansion
  void reset_coefs(void);

  //! Make coefficients after accumulation
  void make_coefs(void);

  //! Accumulate new coefficients
  void accumulate(double x, double y, double z, double mass);
  
  //! Get potential, density and forces at a point
  virtual void getFields
  (double x, double y, double z,
   double& tdens0, double& tpotl0, double& tdens, double& tpotl, 
   double& tpotx, double& tpoty, double& tpotz);


  //! Return current maximum harmonic order in expansion
  int getLmax() { return lmax; }

  //! Return current maximum order in radial expansion
  int getNmax() { return nmax; }
};

/**
   Uses EmpCylSL basis to evaluate expansion coeffients and provide
   potential and density basis fields
*/
class CylindricalSL : public BasicBasis
{

private:

  std::shared_ptr<EmpCylSL> sl;
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

  //! Constructor
  CylindricalSL(const YAML::Node& conf);
  
  //! Destructor
  ~CylindricalSL(void);
  
  //! Prescaling factor
  void set_scale(const double scl) { rscl = scl; }

  //! Zero out coefficients to prepare for a new expansion
  void reset_coefs(void);

  //! Make coefficients after accumulation
  void make_coefs(void);

  //! Accumulate new coefficients
  void accumulate(double x, double y, double z, double mass);
  
  //! Get potential, density and forces at a point
  virtual void getFields
  (double x, double y, double z,
   double& tdens0, double& tpotl0, double& tdens, double& tpotl, 
   double& tpotx, double& tpoty, double& tpotz);


  //! Return current maximum harmonic order in expansion
  int getMmax() { return mmax; }

  //! Return current maximum order in radial expansion
  int getNmax() { return nmax; }
};

//! A factory class for bases
class BasisClient
{
protected:
  //! The wrapped base class instance
  std::shared_ptr<BasicBasis> basis;

public:

  //! Factory constructor
  BasisClient(const YAML::Node& conf);

  //! Get potential, density and forces at a point
  void getFields
  (double x, double y, double z,
   double& tdens0, double& tpotl0, double& tdens, double& tpotl, 
   double& tpotx, double& tpoty, double& tpotz)
  {
    basis->getFields(x, y, z,
		     tdens0, tpotl0, tdens, tpotl,
		     tpotx, tpoty, tpotz);
  }

  //! Accumulate new coefficients
  void accumulate(double x, double y, double z, double mass)
  {
    coefs->accumulate(x, y, z, mass);
  }

  //! Get mass on grid
  double getMass(void) { return basis->getMass(); }

};

#endif // _BasisClient_H


