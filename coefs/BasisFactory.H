#ifndef _BasisFactory_H
#define _BasisFactory_H

#include <functional>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids
#include <yaml-cpp/yaml.h>

#include <ParticleReader.H>
#include <Coefficients.H>
#include <SLGridMP2.H>
#include <YamlCheck.H>
#include <BiorthCyl.H>
#include <EmpCylSL.H>
#include <localmpi.H>

namespace Field
{
  class FieldGenerator;
}

namespace BasisClasses
{
  using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    //! Callback function signature for selection particles to
    //! accumulate
    using Callback =
      std::function<bool(
			 double,	       // mass
			 std::vector<double>&, // position
			 std::vector<double>&, // velocity
			 uint64_t	       // index
			 )>;

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

    //! The callback particle selector
    Callback ftor;

    //! Return instance for coefficients created from arrays
    CoefClasses::CoefStrPtr coefret;

    //! Current particle index
    unsigned long coefindx;
    
  protected:
    
    //! The expansion center
    std::vector<double> coefctr;

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
    
    //! Evaluate fields in spherical coordinates in centered coordinate system
    virtual void all_eval(double r, double costh, double phi,
			  double& den0, double& den1,
			  double& pot0, double& pot1,
			  double& potr, double& pott, double& potp) = 0;
    
    //! Name of the derived class
    std::string name;
    
    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time) = 0;

    //! Using MPI
    bool use_mpi;

    //! Return readable class name
    virtual const std::string classname() = 0;

    //! Subspace index
    virtual const std::string harmonic() = 0;

  public:
    
    //! Constructor from YAML node
    Basis(const YAML::Node& conf);
    
    //! Constructor from YAML string
    Basis(const std::string& confstr);
    
    virtual ~Basis(void) {}
    
    //! Generate coeffients from a particle reader and optional center
    //! location for the expansion
    CoefClasses::CoefStrPtr createFromReader
    (PR::PRptr reader, std::vector<double> center={0.0, 0.0, 0.0});
    
    //! Generate coeffients from an array and optional center location
    //! for the expansion
    CoefClasses::CoefStrPtr createFromArray
    (Eigen::VectorXd& m, RowMatrixXd& p, double time=0.0,
     std::vector<double> center={0.0, 0.0, 0.0}, bool roundrobin=true);
    
    //! Generate coeffients from an array and optional center location
    //! for the expansion using multiple particle partitions
    //@{

    //! Initialize accumulating coefficients from arrays with an optional
    //! center vector.  This is called once to initialize the accumulation.
    void initFromArray
    (std::vector<double> center={0.0, 0.0, 0.0});
    
    //! Initialize accumulating coefficients from arrays.  This is
    //! called once to initialize the accumulation.
    void addFromArray
    (Eigen::VectorXd& m, RowMatrixXd& p, bool roundrobin=true);
    
    //! Create and the coefficients from the array accumulation with the
    //! provided time value
    CoefClasses::CoefStrPtr makeFromArray(double time=0.0);
    //@}
    
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
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs) = 0;

    //@{
    //! The coefficient factory.  Creates a shared pointer to derived
    //! instance of the desired type and returns it.
    static std::shared_ptr<Basis> factory(const YAML::Node&  conf);
    static std::shared_ptr<Basis> factory_string(const std::string& conf);
    //@}
    
    //! Set the particle selector callback
    void setSelector(const Callback& f) { ftor = f; }

    //! Clear the particle selector callback
    void clrSelector() { ftor = nullptr; }
  };
  
  /**
     Uses SLGridSph basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class SphericalSL : public Basis
  {
    
  public:

    using BasisMap   = std::map<std::string, Eigen::VectorXd>;
    using BasisArray = std::vector<std::vector<BasisMap>>;

  private:
    
    //! Helper for constructor
    void initialize();

    std::shared_ptr<SLGridSph> sl;
    std::shared_ptr<SphericalModelTable> mod;
    
    std::string model_file;
    int lmax, nmax, cmap, numr;
    double rmin, rmax, rmap;
    
    bool NO_L0, NO_L1, EVEN_L, EVEN_M, M0_only;
    
    std::vector<Eigen::MatrixXd> potd, dpot, dpt2, dend;
    std::vector<Eigen::MatrixXd> legs, dlegs, d2legs;
    
    Eigen::MatrixXd factorial;
    Eigen::MatrixXd expcoef;
    double scale;
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
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "SphericalSL";}

    //! Subspace index
    virtual const std::string harmonic() { return "l";}

  public:
    
    //! Constructor from YAML node
    SphericalSL(const YAML::Node& conf);
    
    //! Constructor from YAML string
    SphericalSL(const std::string& confstr);
    
    //! Destructor
    virtual ~SphericalSL(void) {}
    
    //! Print and return the cache parameters
    static std::map<std::string, std::string>
    cacheInfo(const std::string& cachefile)
    {
      return SLGridSph::cacheInfo(cachefile);
    }

    //! Prescaling factor
    void set_scale(const double scl) { scale = scl; }
    
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

    //! Return potential-density pair of a vector of a vector of 1d
    //! basis-function grids for SphericalSL, logarithmically spaced
    //! between [logxmin, logxmax] (base 10).
    BasisArray getBasis
    (double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=40)
    {
      return sl->orthoCheck(knots);
    }

  };
  
  /**
     Uses the BiorthCyl basis to evaluate expansion coeffients and
     provide potential and density basis fields
  */
  class FlatDisk : public Basis
  {

  public:

    using BasisMap   = std::map<std::string, Eigen::VectorXd>;
    using BasisArray = std::vector<std::vector<BasisMap>>;

  private:
    
    //! Helper for constructor
    void initialize();

    std::shared_ptr<BiorthCyl> ortho;
    
    std::string model_file;
    int mmax, nmax, cmap, numr;
    double rcylmin, rcylmax;
    
    bool NO_M0, NO_M1, EVEN_M, M0_only;
    
    std::vector<Eigen::MatrixXd> potd, potR, potZ, dend;
    
    Eigen::MatrixXd expcoef;
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

    //! Evaluate basis in cylindrical coordinates
    void all_eval_cyl
    (double R, double z, double phi, 
     double& den0, double& den1,
     double& pot0, double& pot1,
     double& potR, double& potz, double& potp);
      
    //! Evaluate basis in spherical coordinates.  Conversion from the
    //! cylindrical evaluation above.
    virtual void all_eval(double r, double costh, double phi,
			  double& den0, double& den1,
			  double& pot0, double& pot1,
			  double& potr, double& pott, double& potp);
    
    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "FlatDisk";}

    //! Subspace index
    virtual const std::string harmonic() { return "m";}

  public:
    
    //! Constructor from YAML node
    FlatDisk(const YAML::Node& conf);
    
    //! Constructor from YAML string
    FlatDisk(const std::string& confstr);
    
    //! Destructor
    virtual ~FlatDisk(void) {}
    
    //! Print and return the cache parameters
    static std::map<std::string, std::string>
    cacheInfo(const std::string& cachefile)
    {
      return BiorthCyl::cacheInfo(cachefile);
    }

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

    //! Return a vector of a vector of 1d basis-function grids for
    //! FlatDisk, logarithmically spaced between [logxmin, logxmax]
    //! (base 10).
    BasisArray getBasis
    (double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck();

  };
  
  /**
     Uses EmpCylSL basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class Cylindrical : public Basis
  {
    
  public:

    using BasisMap   = std::map<std::string, Eigen::MatrixXd>;
    using BasisArray = std::vector<std::vector<BasisMap>>;

  private:
    
    //! Help for constructor
    void initialize();

    std::shared_ptr<EmpCylSL> sl;
    int lmaxfid, nmaxfid, mmax, mlim, nmax;
    int ncylodd, ncylnx, ncylny, ncylr, cmap, cmapR, cmapZ, vflag;
    int rnum, pnum, tnum;
    double rmin, rmax, rcylmin, rcylmax;
    double acyl, hcyl, hexp;
    bool expcond, logarithmic, density, EVEN_M;
    
    std::vector<Eigen::MatrixXd> potd, dpot, dpt2, dend;
    std::vector<Eigen::MatrixXd> legs, dlegs, d2legs;
    
    bool coefs_defined;
    Eigen::MatrixXd factorial;
    Eigen::MatrixXd expcoef;
    double rscl;
    int used;
    
    std::string eof_file;
    
    using matT = std::vector<Eigen::MatrixXd>;
    using vecT = std::vector<Eigen::VectorXd>;
    
    Eigen::VectorXd work;
    
    //! For coefficient writing
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
    EigenColMajor;
    

    //@{
    //! Basis construction parameters
    double aratio, hratio, dweight, rwidth, ashift, rfactor, rtrunc, ppow;
    bool Ignore, deproject;

    //! DiskType support
    //
    enum DiskType { constant, gaussian, mn, exponential, doubleexpon };
    DiskType DTYPE;
    std::string mtype, dtype, dmodel;

    static const std::map<std::string, DiskType> dtlookup;
    double DiskDens(double R, double z, double phi);
    double dcond(double R, double z, double phi, int M);

  protected:

    //! Evaluate basis in spherical coordinates
    virtual void all_eval(double r, double costh, double phi,
			  double& den0, double& den1,
			  double& pot0, double& pot1,
			  double& potr, double& pott, double& potp);

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "Cylindrical";}

    //! Subspace index
    virtual const std::string harmonic() { return "m";}

  public:
    
    //! Constructor from YAML node
    Cylindrical(const YAML::Node& conf);
    
    //! Constructor from YAML string
    Cylindrical(const std::string& confstr);
    
    //! Destructor
    virtual ~Cylindrical(void) {}
    
    //! Print and return the cache parameters
    static std::map<std::string, std::string>
    cacheInfo(const std::string& cachefile)
    {
      return EmpCylSL::cacheInfo(cachefile);
    }

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
    BasisArray getBasis
    (double xmin=0.0, double xmax=1.0, int numr=40,
     double zmin=-0.1, double zmax=0.1, int numz=40);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck()
    {
      return sl->orthoCheck();
    }
  };
  
  //! Time-dependent potential-density model
  using BasisCoef = std::tuple<std::shared_ptr<Basis>, std::shared_ptr<CoefClasses::Coefs>>;

  //! The function signature
  using AccelFunctor =
    std::function<Eigen::MatrixXd&(double, Eigen::MatrixXd&, Eigen::MatrixXd&, BasisCoef)>;
  
  //! Acceleration function class
  class AccelFunc
  {
  protected:

    //! This evaluates the accelerations for all particles from the
    //! current component.  This is always the same and implemented in
    //! this base class.
    Eigen::MatrixXd&
    evalaccel(Eigen::MatrixXd& ps, Eigen::MatrixXd& accel, BasisCoef mod);

  public:
    
    //! Constructor
    AccelFunc() {}

    //! Destructor
    virtual ~AccelFunc() {}

    //! This evaluates the coefficients at time t.  This is the member
    //! that is implemented in derived classes.  It is more naturally
    //! a protected member in C++, but that makes deriving a class in
    //! Python more difficult.
    virtual void
    evalcoefs(double t, BasisCoef mod) = 0;

    //! This computes and return the acceleration array
    Eigen::MatrixXd&
    F(double t, Eigen::MatrixXd& ps, Eigen::MatrixXd& accel,  BasisCoef mod)
    {
      // Evaluate the coefficients for time t
      evalcoefs(t, mod);

      // Compute and return the accumulated acceleration
      return evalaccel(ps, accel, mod);
    }
  };

  //! Evaluate acceleration for one component, return acceleration
  class AllTimeAccel : public AccelFunc
  {
  public:
    
    //! Constructor
    AllTimeAccel() {}

    //! Interpolate and install coefficients at time t
    void evalcoefs(double t, BasisCoef mod);

  };

  //! Evaluate acceleration for one component, return acceleration
  class SingleTimeAccel : public AccelFunc
  {
  public:
    
    //! Initialize and install coefficients for the given time
    SingleTimeAccel(double time, std::vector<BasisCoef> mod);

    //! No coefficient evaluation
    void evalcoefs(double t, BasisCoef mod) {}

  };


  std::tuple<Eigen::VectorXd, Eigen::Tensor<float, 3>>
  IntegrateOrbits (double tinit, double tfinal, double h,
		   Eigen::MatrixXd ps, std::vector<BasisCoef> bfe,
		   AccelFunctor F, int nout=std::numeric_limits<int>::max());

  using BasisPtr = std::shared_ptr<Basis>;
}
// END: namespace BasisClasses
  
#endif // _BasisFactory_H
  
  
  
