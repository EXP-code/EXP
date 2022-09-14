#ifndef _BasisFactory_H
#define _BasisFactory_H

#include <Eigen/Eigen>
#include <yaml-cpp/yaml.h>

#include <ParticleReader.H>
#include <CoefContainer.H>
#include <SLGridMP2.H>
#include <EmpCylSL.H>
#include <localmpi.H>

namespace Field
{
  class FieldGenerator;
}

namespace Basis
{
    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  /** 
      Abstract class for any type of basis for evaluating density,
      potential and acceleration.
      
      Provides a factory member that returns derived class instances.
  */
  //! Defines a basis-based 
  class Basis
  {
    friend class Field::FieldGenerator;

  private:

    //! Helper for constructor
    void initialize();

    //! Helper for factory
    static std::shared_ptr<Basis> factory_initialize(const YAML::Node& node);

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
    
    //! Evaluate fields in spherical coordinates
    virtual void all_eval(double r, double costh, double phi,
			  double& den0, double& den1,
			  double& pot0, double& pot1,
			  double& potr, double& pott, double& potp) = 0;
    
    //! Name of the derived class
    std::string name;
    
    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(Coefs::CoefStrPtr coefs, double time) = 0;

    //! Using MPI
    bool use_mpi;

  public:
    
    //! Constructor from YAML node
    Basis(const YAML::Node& conf);
    
    //! Constructor from YAML string
    Basis(const std::string& confstr);
    
    virtual ~Basis(void) {}
    
    //! Generate coeffients from a particle reader and optional center
    //! location for the expansion
    Coefs::CoefStrPtr createFromReader
    (PR::PRptr reader, std::vector<double> center={0.0, 0.0, 0.0});
    
    //! Generate coeffients from an array and optional center location
    //! for the expansion
    Coefs::CoefStrPtr createFromArray
    (RowMatrixXd& p, double time=0.0,
     std::vector<double> center={0.0, 0.0, 0.0});
    
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
    
    //! Provide a set of coefficients using a CoefStruct
    virtual void set_coefs(Coefs::CoefStrPtr coefs) = 0;

    //@{
    //! The coefficient factory.  Creates a shared pointer to derived
    //! instance of the desired type and returns it.
    static std::shared_ptr<Basis> factory(const YAML::Node&  conf);
    static std::shared_ptr<Basis> factory_string(const std::string& conf);
    //@}
    

  };
  
  /**
     Uses SLGridSph basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class SphericalSL : public Basis
  {
    
  private:
    
    //! Helper for constructor
    void initialize();

    std::shared_ptr<SLGridSph> sl;
    std::shared_ptr<SphericalModelTable> mod;
    
    std::string model_file;
    int lmax, nmax, cmap, numr;
    double rmin, rmax;
    
    bool NO_L0, NO_L1, EVEN_L, EVEN_M, M0_only;
    
    Eigen::MatrixXd potd, dpot, dpt2, dend;
    Eigen::MatrixXd legs, dlegs, d2legs;
    
    Eigen::MatrixXd factorial;
    Eigen::MatrixXd expcoef;
    double rscl;
    int N1, N2;
    int used;
    
    using matT = std::vector<Eigen::MatrixXd>;
    using vecT = std::vector<Eigen::VectorXd>;
    
    double totalMass;
    int npart;
    
    Eigen::VectorXd work;
    
    //! For coefficient writing
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
    EigenColMajor;
    
  protected:

    //! Evaluate basis in spherical coordinates
    virtual void all_eval(double r, double costh, double phi,
			  double& den0, double& den1,
			  double& pot0, double& pot1,
			  double& potr, double& pott, double& potp);
    
    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(Coefs::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(Coefs::CoefStrPtr coefs);

  public:
    
    //! Constructor from YAML node
    SphericalSL(const YAML::Node& conf);
    
    //! Constructor from YAML string
    SphericalSL(const std::string& confstr);
    
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

    //! Return a vector of a vector of 1d basis-function grids for
    //! SphericalSL, logarithmically spaced between [logxmin, logxmax]
    //! (base 10).
    std::vector<std::vector<Eigen::VectorXd>> getBasis
    (double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);
  };
  
  /**
     Uses EmpCylSL basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class Cylindrical : public Basis
  {
    
  private:
    
    //! Help for constructor
    void initialize();

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

    //! Evaluate basis in spherical coordinates
    virtual void all_eval(double r, double costh, double phi,
			  double& den0, double& den1,
			  double& pot0, double& pot1,
			  double& potr, double& pott, double& potp);

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(Coefs::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(Coefs::CoefStrPtr coefs);

  public:
    
    //! Constructor from YAML node
    Cylindrical(const YAML::Node& conf);
    
    //! Constructor from YAML string
    Cylindrical(const std::string& confstr);
    
    //! Destructor
    virtual ~Cylindrical(void) {}
    
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

    //! Return a vector of vector of basis functions for Cylindrical,
    //! linearly spaced between [xmin, xmax] and [zmin, zmax]
    std::vector<std::vector<Eigen::MatrixXd>> getBasis
    (double xmin=0.0, double xmax=1.0, int numr=40,
     double zmin=-0.1, double zmax=0.1, int numz=40);
  };
  
  using BasisPtr = std::shared_ptr<Basis>;
}
  
#endif // _BasisFactory_H
  
  
  
