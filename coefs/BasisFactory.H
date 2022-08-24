#ifndef _BasisFactory_H
#define _BasisFactory_H

#include <Eigen/Eigen>
#include <yaml-cpp/yaml.h>

#include <ParticleReader.H>
#include <CoefFactory.H>
#include <SLGridMP2.H>
#include <EmpCylSL.H>
#include <localmpi.H>

namespace Basis
{
  
  //! Defines a basis-based potential and acceleration class
  class Basis
  {
  protected:
    
    //! Contains contructor and BFE parameter database
    YAML::Node node, conf;
    
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
    
    void all_eval(double r, double costh, double phi,
		  double& den0, double& den1,
		  double& pot0, double& pot1,
		  double& potr, double& pott, double& potp);
    
    //! Name of the derived class
    std::string name;
    
    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(Coefs::CoefStrPtr coefs, double time) = 0;

  public:
    
    //! Using MPI
    static bool mpi;
    
    //! Constructor
    Basis(const YAML::Node& conf);
    
    virtual ~Basis(void) {}
    
    //! Generate coeffients from a particle reader
    Coefs::CoefStrPtr createCoefficients(PR::PRptr reader);
    
    //! Get potential, density and forces at a point
    virtual void getFields
    (double x, double y, double z,
     double& tdens0, double& tpotl0, double& tdens, double& tpotl, 
     double& tpotx, double& tpoty, double& tpotz) = 0;
    
    //! Accumulate new coefficients
    virtual void accumulate(double x, double y, double z, double mass) = 0;
    
    //! Get mass on grid
    double getMass(void) { return totalMass; }
    
    //! Zero out coefficients to prepare for a new expansion
    virtual void reset_coefs(void) = 0;
    
    //! Make coefficients after accumulation
    virtual void make_coefs(void) = 0;
    
    //! The coefficient factory.  Creates a shared pointer to derived
    //! instance of the desired type and returns it.
    static std::shared_ptr<Basis> factory(const YAML::Node& conf);
    
  };
  
  /**
     Uses SLGridSph basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class SphericalSL : public Basis
  {
    
  private:
    
    std::shared_ptr<SLGridSph> sl;
    std::shared_ptr<SphericalModelTable> mod;
    
    std::string model_file;
    int lmax, nmax, cmap, numr;
    double rmin, rmax;
    
    bool NO_L0, NO_L1, EVEN_L, EVEN_M, M0_only;
    
    Eigen::MatrixXd potd, dpot, dpt2, dend;
    Eigen::MatrixXd legs, dlegs, d2legs;
    
    bool coefs_defined;
    Eigen::MatrixXd factorial;
    Eigen::MatrixXd expcoef;
    double rscl;
    int N1, N2;
    int used;
    
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
    
    //! Evaluate basis in spherical coordinates
    void all_eval(double r, double costh, double phi,
		  double& den0, double& den1,
		  double& pot0, double& pot1,
		  double& potr, double& pott, double& potp);
    
    //! For coefficient writing
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
    EigenColMajor;
    
  protected:

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(Coefs::CoefStrPtr coefs, double time);

  public:
    
    //! Constructor
    SphericalSL(const YAML::Node& conf);
    
    //! Destructor
    virtual ~SphericalSL(void) {}
    
    //! Prescaling factor
    void set_scale(const double scl) { rscl = scl; }
    
    //! Zero out coefficients to prepare for a new expansion
    void reset_coefs(void);
    
    //! Make coefficients after accumulation
    void make_coefs(void);
    
    //! Accumulate new coefficients
    virtual void accumulate(double x, double y, double z, double mass);
    
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
  class CylindricalSL : public Basis
  {
    
  private:
    
    std::shared_ptr<EmpCylSL> sl;
    int lmax, mmax, nmax, mlim;
    int ncylorder, ncylodd, ncylnx, ncylny, ncylr, cmap, cmapR, cmapZ, vflag;
    int rnum, pnum, tnum, ashift;
    double rmin, rmax, rcylmin, rcylmax;
    double acyl, hcyl, hexp;
    bool expcond, logarithmic, density, EVEN_M;
    
    Eigen::MatrixXd potd, dpot, dpt2, dend;
    Eigen::MatrixXd legs, dlegs, d2legs;
    
    bool coefs_defined;
    Eigen::MatrixXd factorial;
    Eigen::MatrixXd expcoef;
    double rscl;
    int used;
    
    bool eof_over, try_cache;
    std::string eof_file;
    
    using matT = std::vector<Eigen::MatrixXd>;
    using vecT = std::vector<Eigen::VectorXd>;
    
    Eigen::VectorXd work;
    
    //! For coefficient writing
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
    EigenColMajor;
    
  protected:

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(Coefs::CoefStrPtr coefs, double time);

  public:
    
    //! Constructor
    CylindricalSL(const YAML::Node& conf);
    
    //! Destructor
    virtual ~CylindricalSL(void) {}
    
    //! Prescaling factor
    void set_scale(const double scl) { rscl = scl; }
    
    //! Zero out coefficients to prepare for a new expansion
    void reset_coefs(void);
    
    //! Make coefficients after accumulation
    void make_coefs(void);
    
    //! Accumulate new coefficients
    virtual void accumulate(double x, double y, double z, double mass);
    
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
  
  using BasisPtr = std::shared_ptr<Basis>;
}
  
#endif // _BasisFactory_H
  
  
  
