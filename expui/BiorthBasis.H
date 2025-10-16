#ifndef _BiorthBasis_H
#define _BiorthBasis_H

#include <functional>
#include <tuple>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids
#include <yaml-cpp/yaml.h>

#include <DiskDensityFunc.H>
#include <ParticleReader.H>
#include <Coefficients.H>
#include <BiorthBess.H>
#include <BasisFactory.H>
#include <BiorthCube.H>
#include <sltableMP2.H>
#include <SLGridMP2.H>
#include <YamlCheck.H>
#include <BiorthCyl.H>
#include <EmpCylSL.H>
#include <localmpi.H>
#include <exputils.H>

namespace Field
{
  class FieldGenerator;
}

namespace BasisClasses
{
  /** 
      Abstract class for any type of basis for evaluating density,
      potential and acceleration.
      
      Provides a factory member that returns derived class instances.
  */
  //! Defines a basis-based 
  class BiorthBasis : public Basis
  {
    friend class Field::FieldGenerator;

  private:

    //! Helper for factory
    static std::shared_ptr<BiorthBasis>
    factory_initialize(const YAML::Node& node);

  protected:
    
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
    
    //! Get field labels
    std::vector<std::string> getFieldLabels(const Coord ctype);

    //! Evaluate fields in spherical coordinates in centered coordinate system
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi) = 0;

    //! Evaluate fields in cylindrical coordinates in centered coordinate system
    virtual std::vector<double>
    cyl_eval(double r, double costh, double phi) = 0;


    //! Evaluate fields in Cartesian coordinates in centered coordinate system
    virtual std::vector<double>
    crt_eval(double x, double y, double z) = 0;

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time) = 0;

    //! Subspace index
    virtual const std::string harmonic() = 0;

    //@{
    //! Compute acceleration helper for each derived class
    virtual void computeAccel(double x, double y, double z,
			      Eigen::Ref<Eigen::Vector3d> vstore) = 0;

    //! C++ local storage for computeAccel
    Eigen::Vector3d vstore;

    //! C++ local storage for getAccelArray
    RowMatrixXd varray;
    //@}

  public:
    
    //! Constructor from YAML node
    BiorthBasis(const YAML::Node& conf,
		const std::string& name="BiorthBasis") : Basis(conf, name) {}
    
    //! Constructor from YAML string
    BiorthBasis(const std::string& confstr,
		const std::string& name="BiorthBasis") : Basis(confstr, name) {}
    
    //! Destructor
    virtual ~BiorthBasis(void) {}
    
    //! Generate coeffients from a particle reader and optional center
    //! location for the expansion
    CoefClasses::CoefStrPtr createFromReader
    (PR::PRptr reader, Eigen::Vector3d center=Eigen::Vector3d::Zero(),
     RowMatrix3d rot=RowMatrix3d::Identity());
    
    //! Generate coeffients from an array and optional center location
    //! for the expansion
    CoefClasses::CoefStrPtr createFromArray
    (Eigen::VectorXd& m, RowMatrixXd& p, double time=0.0,
     Eigen::Vector3d center=Eigen::Vector3d::Zero(),
     RowMatrix3d rot=RowMatrix3d::Identity(),
     bool roundrobin=true, bool posvelrows=false);
    
    //! Generate coeffients from an array and optional center location
    //! for the expansion using multiple particle partitions
    //@{

    //! Initialize accumulating coefficients from arrays with an optional
    //! center vector.  This is called once to initialize the accumulation.
    void initFromArray
    (Eigen::Vector3d center=Eigen::Vector3d::Zero(),
     RowMatrix3d rot=RowMatrix3d::Identity());
    
    //! Initialize accumulating coefficients from arrays.  This is
    //! called once to initialize the accumulation.
    void addFromArray
    (Eigen::VectorXd& m, RowMatrixXd& p,
     bool roundrobin=true, bool posvelrows=false);
    
    //! Create and the coefficients from the array accumulation with the
    //! provided time value
    CoefClasses::CoefStrPtr makeFromArray(double time=0.0);
    //@}
    
    //! Create and the coefficients from a function callback with the
    //! provided time value.  Set potential=true if function is a
    //! potential field.  Density is assumed if false (default).
    virtual CoefClasses::CoefStrPtr makeFromFunction
    (std::function<double(double x, double y, double z, double time)> func,
     std::map<std::string, double>& params, double time=0.0, bool potential=false)
    {
      throw std::runtime_error("Error: BiorthBasis::makeFromFunction is not implemented");
      return nullptr;
    }
    
    //! Compute a quadrature over the basis domain using function
    //! callback for the integrand
    virtual double computeQuadrature
    (std::function<double(double x, double y, double z)> func,
     std::map<std::string, double>& params)
    {
      throw std::runtime_error("Error: BiorthBasis::computeQuadrature is not implemented");
      return 0.0;
    }
    
    //! Accumulate new coefficients
    virtual void accumulate(double x, double y, double z, double mass) = 0;
    
    //! Accumulate new coefficients
    virtual void accumulate(double mass,
			    double x, double y, double z,
			    double u, double v, double w)
    { accumulate(x, y, z, mass); }
    
    //! Get mass on grid
    double getMass(void) { return totalMass; }
    
    //! Zero out coefficients to prepare for a new expansion
    virtual void reset_coefs(void) = 0;
    
    //! Make coefficients after accumulation
    virtual void make_coefs(void) = 0;
    
    //! Provide a set of coefficients using a CoefStruct
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs) = 0;

    //! Set field coordindate system
    void setFieldType(std::string coord_type)
    { coordinates = parseFieldType(coord_type); }
    
    //! Get current field coordinate type
    std::string getFieldType() { return coordLabels[coordinates]; }
    
    //@{
    //! The coefficient factory.  Creates a shared pointer to derived
    //! instance of the desired type and returns it.
    static std::shared_ptr<BiorthBasis> factory(const YAML::Node&  conf);
    static std::shared_ptr<BiorthBasis> factory_string(const std::string& conf);
    //@}
    
    //! Set the particle selector callback
    void setSelector(const Callback& f) { ftor = f; }

    //! Clear the particle selector callback
    void clrSelector() { ftor = nullptr; }

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system
    virtual Eigen::Vector3d& getAccel(double x, double y, double z)
    {
      computeAccel(x, y, z, vstore);
      return vstore;
    }

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system for a collections of points
    virtual RowMatrixXd& getAccelArray(Eigen::VectorXd& x,
				       Eigen::VectorXd& y,
				       Eigen::VectorXd& z)
    {
      // Sanity check
      if (x.size() != y.size() || x.size() != z.size())
	throw std::runtime_error("BiorthBasis::getAccelArray: "
				 "x, y, z vectors must be of the same size");
      
      // Resize the return array
      varray.resize(x.size(), 3);
#pragma omp parallel for
      for (int i=0; i<x.size(); ++i) {
	computeAccel(x(i), y(i), z(i), varray.row(i));
      }

      return varray;
    }
  };

  /**
     An abstract spherical basis to evaluate expansion coeffients and
     provide potential and density basis fields
  */
  class Spherical : public BiorthBasis
  {
    
  public:

    using BasisMap   = std::map<std::string, Eigen::VectorXd>;
    using BasisArray = std::vector<std::vector<BasisMap>>;

  protected:
    
    //! Helper for constructor
    void initialize();

    //! Load coefficients for a particular time
    virtual void load_coefs(CoefClasses::CoefStrPtr coef, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "Spherical";}

    //! Subspace index
    virtual const std::string harmonic() { return "l";}

    //! Evaluate basis in spherical coordinates
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi);

    //! Evaluate basis in cartesian coordinates
    virtual std::vector<double>
    crt_eval(double x, double y, double z);

    //! Evaluate basis in cylindrical coordinates
    virtual std::vector<double>
    cyl_eval(double R, double z, double phi);

    //@{
    //! Required basis members

    //! Get potential
    virtual void get_pot(Eigen::MatrixXd& tab, double x) = 0;

    //! Get density
    virtual void get_dens(Eigen::MatrixXd& tab, double x) = 0;

    //! Get force
    virtual void get_force(Eigen::MatrixXd& tab, double x) = 0;

    //@}

    //@{
    //! Internal parameters and storage
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

    //@}

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system
    virtual void computeAccel(double x, double y, double z, 
			      Eigen::Ref<Eigen::Vector3d> vstore);

  public:
    
    //! Constructor from YAML node
    Spherical(const YAML::Node& conf, const std::string& forceID);
    
    //! Constructor from YAML string
    Spherical(const std::string& confstr, const std::string& forceID);
    
    //! Destructor
    virtual ~Spherical(void) {}
    
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
    
    //! Return current maximum harmonic order in expansion
    int getLmax() { return lmax; }
    
    //! Return current maximum order in radial expansion
    int getNmax() { return nmax; }

    //! Return potential-density pair of a vector of a vector of 1d
    //! basis-function grids for Spherical [rmin, rmax]
    virtual BasisArray getBasis
    (double rmin=0.0, double rax=1.0, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    virtual std::vector<Eigen::MatrixXd> orthoCheck(int knots) = 0;

    //! Biorthogonality sanity check
    bool orthoTest(int knots=100)
    {
      auto [ret, worst, lworst] = orthoCompute(orthoCheck(knots));
      // For the CTest log
      if (myid==0)
	std::cout << "---- Spherical::orthoTest: worst=" << worst << std::endl;
      return ret;
    }

    //! Create and the coefficients from a function callback with the
    //! provided time value.  Set potential=true if function is a
    //! potential field.  Density is assumed if false (default).
    virtual CoefClasses::CoefStrPtr makeFromFunction
    (std::function<double(double x, double y, double z, double time)> func,
     std::map<std::string, double>& params, double time=0.0, bool potential=false);


    //! Compute a quadrature over the basis domain using function
    //! callback for the integrand
    double computeQuadrature
    (std::function<double(double x, double y, double z)> func,
     std::map<std::string, double>& params);
  };

  /**
     Uses SLGridSph basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class SphericalSL : public Spherical
  {
    
  protected:
    
    //! Helper for constructor
    void initialize();

    static const std::set<std::string> valid_keys;

    std::shared_ptr<SLGridSph> sl;
    std::shared_ptr<SphericalModelTable> mod;
    
    std::string model_file;
    
    //! Return readable class name
    const std::string classname() { return "SphericalSL";}

    // Get potential
    void get_pot(Eigen::MatrixXd& tab, double r)
    { sl->get_pot(tab, r); }

    // Get density
    void get_dens(Eigen::MatrixXd& tab, double r)
    { sl->get_dens(tab, r); }

    // Get force
    void get_force(Eigen::MatrixXd& tab, double r)
    { sl->get_force(tab, r); }

  public:
    
    //! Constructor from YAML node
    SphericalSL(const YAML::Node& conf);
    
    //! Constructor from YAML string
    SphericalSL(const std::string& confstr);
    
    //! Destructor
    virtual ~SphericalSL(void) {}
    
    //! Return potential-density pair of a vector of a vector of 1d
    //! basis-function grids for SphericalSL, logarithmically spaced
    //! between [logxmin, logxmax] (base 10).
    BasisArray getBasis(double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=40)
    { return sl->orthoCheck(knots); }

  };
  
  /**
     Uses Bessel basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class Bessel : public Spherical
  {
    
  protected:
    
    //! Helper for constructor
    void initialize();

    static const std::set<std::string> valid_keys;

    //! Return readable class name
    const std::string classname() { return "Bessel";}

    //! Grid size for Bessel function table
    int rnum = 2000;

    //! Biorthgonal Bessel function generator
    std::shared_ptr<BiorthBess> bess;
    
    // Get potential
    void get_pot(Eigen::MatrixXd& tab, double r)
    { bess->get_potl(r, tab); }

    // Get density
    void get_dens(Eigen::MatrixXd& tab, double r)
    { bess->get_dens(r, tab); }

    // Get force
    void get_force(Eigen::MatrixXd& tab, double r)
    { bess->get_dpotl(r, tab); }

  public:
    
    //! Constructor from YAML node
    Bessel(const YAML::Node& conf);
    
    //! Constructor from YAML string
    Bessel(const std::string& confstr);
    
    //! Destructor
    virtual ~Bessel(void) {}

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=40)
    { return bess->orthoCheck(knots); }
  };
  

  /**
     Uses the BiorthCyl basis to evaluate expansion coeffients and
     provide potential and density basis fields
  */
  class FlatDisk : public BiorthBasis
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
    virtual std::vector<double>
    cyl_eval(double R, double z, double phi);
      
    //! Evaluate basis in spherical coordinates.  Conversion from the
    //! cylindrical evaluation above.
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi);
    
    // Cartesian
    virtual std::vector<double>
    crt_eval(double x, double y, double z);
    
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

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system
    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> vstore);

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
    
    //! Return current maximum harmonic order in expansion
    int getMmax() { return mmax; }
    
    //! Return current maximum order in radial expansion
    int getNmax() { return nmax; }

    //! Return a vector of a vector of 1d basis-function grids for
    //! FlatDisk, logarithmically spaced between [logxmin, logxmax]
    //! (base 10).
    BasisArray getBasis(double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=400);

    //! Biorthogonality sanity check
    bool orthoTest()
    {
      auto [ret, worst, lworst] = orthoCompute(orthoCheck());
      // For the CTest log
      if (myid==0)
	std::cout << "---- FlatDisk::orthoTest: worst=" << worst << std::endl;
      return ret;
    }

  };
  
  /**
     Uses the Clutton-Brock basis to evaluate expansion coeffients and
     provide potential and density basis fields
  */
  class CBDisk : public BiorthBasis
  {

  public:

    using BasisMap   = std::map<std::string, Eigen::VectorXd>;
    using BasisArray = std::vector<std::vector<BasisMap>>;

  private:
    
    //! Helper for constructor
    void initialize();

    int mmax, nmax;
    double scale;
    
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
    
    //! Get potential
    void get_pot(Eigen::MatrixXd& tab, double r);

    //! Get density
    void get_dens(Eigen::MatrixXd& tab, double r);

    //! Get force
    void get_force(Eigen::MatrixXd& tab, double r);

    //@{
    //! 2D Clutton-Brock basis
    double phif(const int n, const int m, const double r);
    double dphi(const int n, const int m, const double r);
    
    double potl(const int n, const int m, const double r);
    double dpot(const int n, const int m, const double r);
    double dens(const int n, const int m, const double r);
    
    double norm(const int n, const int m);

    void potl(const int m, const double r, Eigen::VectorXd& a);
    void dpot(const int m, const double r, Eigen::VectorXd& a);
    void dens(const int m, const double r, Eigen::VectorXd& a);
    //@}

    //! Potential, force, and density scaling
    double fac1, fac2;

  protected:

    //! Evaluate basis in cylindrical coordinates
    virtual std::vector<double>
    cyl_eval(double R, double z, double phi);
      
    //! Evaluate basis in spherical coordinates.  Conversion from the
    //! cylindrical evaluation above.
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi);
    
    // Cartesian
    virtual std::vector<double>
    crt_eval(double x, double y, double z);
    
    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "CBDisk";}

    //! Subspace index
    virtual const std::string harmonic() { return "m";}

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system
    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> vstore);

  public:
    
    //! Constructor from YAML node
    CBDisk(const YAML::Node& conf);
    
    //! Constructor from YAML string
    CBDisk(const std::string& confstr);
    
    //! Destructor
    virtual ~CBDisk(void) {}
    
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
    
    //! Return current maximum harmonic order in expansion
    int getMmax() { return mmax; }
    
    //! Return current maximum order in radial expansion
    int getNmax() { return nmax; }

    //! Return a vector of a vector of 1d basis-function grids for
    //! CBDisk, logarithmically spaced between [logxmin, logxmax]
    //! (base 10).
    BasisArray getBasis(double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=4000);

    //! Biorthogonality sanity check
    bool orthoTest()
    {
      auto [ret, worst, lworst] = orthoCompute(orthoCheck());
      // For the CTest log
      if (myid==0)
	std::cout << "CBDisk::orthoTest: worst=" << worst << std::endl;
      return ret;
    }

  };
  
  /**
     Uses EmpCylSL basis to evaluate expansion coeffients and provide
     potential and density basis fields
  */
  class Cylindrical : public BiorthBasis
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
    double acyl, hcyl;
    bool expcond, logarithmic, density, EVEN_M, sech2 = false;
    
    std::vector<Eigen::MatrixXd> potd, dpot, dpt2, dend;
    std::vector<Eigen::MatrixXd> legs, dlegs, d2legs;
    
    bool coefs_defined;
    Eigen::MatrixXd factorial;
    Eigen::MatrixXd expcoef;
    double rscl;
    int used;
    
    std::string cachename;
    bool oldcache = false;
    
    using matT = std::vector<Eigen::MatrixXd>;
    using vecT = std::vector<Eigen::VectorXd>;
    
    Eigen::VectorXd work;
    
    //! For coefficient writing
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
    EigenColMajor;
    

    //@{
    //! Basis construction parameters
    double aratio, hratio, dweight, rwidth, ashift, rfactor, rtrunc, ppow, Mfac, HERNA;
    bool Ignore, deproject;
    std::string pyname;

    //! DiskType support
    //
    enum DiskType { constant, gaussian, mn, exponential, doubleexpon, diskbulge, python };
    DiskType DTYPE;
    std::string mtype, dtype, dmodel;

    static const std::map<std::string, DiskType> dtlookup;
    double DiskDens(double R, double z, double phi);
    double dcond(double R, double z, double phi, int M);
    std::shared_ptr<DiskDensityFunc> pyDens;

  protected:

    //! Evaluate basis in spherical coordinates
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi);

    //! Evaluate basis in cylindrical coordinates
    virtual std::vector<double>
    cyl_eval(double r, double z, double phi);

    //! Evaluate basis in cartesian coordinates
    virtual std::vector<double>
    crt_eval(double x, double y, double z);

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

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system
    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> vstore);

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
    
    //! Return current maximum harmonic order in expansion
    int getMmax() { return mmax; }
    
    //! Return current maximum order in radial expansion
    int getNmax() { return nmax; }

    //! Return a vector of vector of basis functions for Cylindrical,
    //! linearly spaced between [xmin, xmax] and [zmin, zmax]
    BasisArray getBasis
    (double xmin=0.0, double xmax=1.0, int numr=40,
     double zmin=-0.1, double zmax=0.1, int numz=40, bool linear=true);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=40)
    {
      return sl->orthoCheck();
    }

    //! Biorthogonality sanity check
    bool orthoTest()
    {
      auto [ret, worst, lworst] = orthoCompute(sl->orthoCheck());
      // For the CTest log
      if (myid==0)
	std::cout << "---- Cylindrical::orthoTest: worst=" << worst << std::endl;
      return ret;
    }

    //! Create and the coefficients from a function callback with the
    //! provided time value.  Set potential=true if function is a
    //! potential field.  Density is assumed if false (default).
    virtual CoefClasses::CoefStrPtr makeFromFunction
    (std::function<double(double x, double y, double z, double time)> func,
     std::map<std::string, double>& params, double time=0.0, bool potential=false);

    //! Compute a quadrature over the basis domain using function
    //! callback for the integrand
    double computeQuadrature
    (std::function<double(double x, double y, double z)> func,
     std::map<std::string, double>& params);
  };
  
  /**
     Uses the SLGridSlab basis to evaluate expansion coeffients and
     provide potential and density basis fields
  */
  class Slab : public BiorthBasis
  {

  public:
    
    using BasisMap   = std::map<std::string, Eigen::VectorXd>;
    using BasisArray = std::vector<std::vector<std::vector<BasisMap>>>;

  private:
    
    //! Wave function constant
    static constexpr std::complex<double> kfac{0.0, 2.0*M_PI};

    //! Initialization helper
    void initialize();

    //! Orthogonal basis instance
    std::shared_ptr<SLGridSlab> ortho;
    
    //! Minimum expansion order for restriction
    int nminx, nminy, nminz;

    //! Maximum expansion order for construction
    int nmaxx, nmaxy, nmaxz, imx, imy, imz;

    //! Number of integration knots for orthogonal check
    int knots;

    //! SLGridSlab mesh size
    int ngrid = 1000;

    //! Target model type for SLGridSlab
    std::string type = "isothermal";

    //! Scale height for slab
    double hslab;

    //! Upper vertical bound for slab
    double zmax;

    //! Number of particles used to compute grid
    unsigned used;
  
    using coefType = Eigen::Tensor<std::complex<double>, 3>;

    coefType expcoef;
    
    //! Notal mass on grid
    double totalMass;

    //! Number of particles
    int npart;
    
  protected:

    //! Evaluate basis in Cartesian coordinates
    virtual std::vector<double>
    crt_eval(double x, double y, double z);
    
    //! Evaluate basis in spherical coordinates.  Conversion from the
    //! Cartesian evaluation above.
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi);
    
    //! Evaluate basis in cylindrical coordinates
    virtual std::vector<double>
    cyl_eval(double r, double z, double phi);

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "Slab";}

    //! Readable index name
    virtual const std::string harmonic()  { return "n";}

    //! Evaluate field
    std::tuple<double, double, double, double, double>
    eval(double x, double y, double z);

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system
    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> vstore);

  public:
    
    //! Constructor from YAML node
    Slab(const YAML::Node& conf);
    
    //! Constructor from YAML string
    Slab(const std::string& confstr);
    
    //! Destructor
    virtual ~Slab(void) {}
    
    //! Zero out coefficients to prepare for a new expansion
    void reset_coefs(void);
    
    //! Make coefficients after accumulation
    void make_coefs(void);
    
    //! Accumulate new coefficients
    virtual void accumulate(double x, double y, double z, double mass);
    
    //! Return current maximum harmonic order in expansion
    Eigen::Vector3i getNmax() { return {nmaxx, nmaxy, nmaxz}; }
    
    //! Return potential-density pair of a vector of a vector of 1d
    //! basis-function grids for the SL slab, linearly spaced
    //! between [zmin, zmax]
    BasisArray getBasis
    (double zmin=-1.0, double zmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=40);

    //! Biorthogonality sanity check
    bool orthoTest()
    {
      auto ret = orthoCheck();
      double worst = 0.0;
      for (auto c : ret) {
	for (int n1=0; n1<c.rows(); n1++) {
	  for (int n2=0; n2<c.cols(); n2++) {
	    double val = std::abs(c(n1, n2));
	    if (n1==n2) val = std::fabs(val - 1.0);
	    worst = std::max<double>(val, worst);
	  }
	}
      }

      if (myid==0)
      std::cout << "---- Slab::orthoTest: worst=" << worst << std::endl;

      if (worst > __EXP__::orthoTol) return false;
      return true;
    }

  };
  
  /**
     Uses the BiorthCyl basis to evaluate expansion coeffients and
     provide potential and density basis fields
  */
  class Cube : public BiorthBasis
  {

  public:
    
    using BasisMap   = std::map<std::string, Eigen::VectorXd>;
    using BasisArray = std::vector<std::vector<BasisMap>>;

  private:
    
    //! Wave function constant
    static constexpr std::complex<double> kfac{0.0, 2.0*M_PI};

    //! Initialization helper
    void initialize();

    //! Orthogonal basis instance
    std::shared_ptr<BiorthCube> ortho;
    
    //! Minimum expansion order for restriction
    int nminx, nminy, nminz;

    //! Maximum expansion order for construction
    int nmaxx, nmaxy, nmaxz;

    //! Number of integration knots for orthogonal check
    int knots;

    //! Number of particles used to compute grid
    unsigned used;
  
    using coefType = Eigen::Tensor<std::complex<double>, 3>;

    coefType expcoef;
    
    //! Notal mass on grid
    double totalMass;

    //! Number of particles
    int npart;
    
  protected:

    //! Evaluate basis in Cartesian coordinates
    virtual std::vector<double>
    crt_eval(double x, double y, double z);
    
    //! Evaluate basis in spherical coordinates.  Conversion from the
    //! Cartesian evaluation above.
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi);
    
    //! Evaluate basis in cylindrical coordinates
    virtual std::vector<double>
    cyl_eval(double r, double z, double phi);

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "Cube";}

    //! Readable index name
    virtual const std::string harmonic()  { return "n";}

    //! Evaluate acceleration in Cartesian coordinates in centered
    //! coordinate system
    void computeAccel(double x, double y, double z, 
		      Eigen::Ref<Eigen::Vector3d> vstore);


  public:
    
    //! Constructor from YAML node
    Cube(const YAML::Node& conf);
    
    //! Constructor from YAML string
    Cube(const std::string& confstr);
    
    //! Destructor
    virtual ~Cube(void) {}
    
    //! Zero out coefficients to prepare for a new expansion
    void reset_coefs(void);
    
    //! Make coefficients after accumulation
    void make_coefs(void);
    
    //! Accumulate new coefficients
    virtual void accumulate(double x, double y, double z, double mass);
    
    //! Return current maximum harmonic order in expansion
    Eigen::Vector3i getNmax() { return {nmaxx, nmaxy, nmaxz}; }
    
    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=40);

    //! Biorthogonality sanity check
    bool orthoTest()
    {
      auto c = orthoCheck()[0];
      double worst = 0.0;
      for (int n1=0; n1<c.rows(); n1++) {
	for (int n2=0; n2<c.cols(); n2++) {
	  double val = std::abs(c(n1, n2));
	  if (n1==n2) val = std::fabs(val - 1.0);
	  worst = std::max<double>(val, worst);
	}
      }
	
      if (myid==0)
	std::cout << "---- Cube::orthoTest: worst=" << worst << std::endl;

      if (worst > __EXP__::orthoTol) return false;
      return true;
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
		   AccelFunctor F, int nout=0);

  using BiorthBasisPtr = std::shared_ptr<BiorthBasis>;
}
// END: namespace BasisClasses
  
#endif // _BiorthBasis_H
  
  
  
