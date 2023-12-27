#ifndef _OrthoBasisFactory_H
#define _OrthoBasisFactory_H

#include <functional>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids
#include <yaml-cpp/yaml.h>

#include <ParticleReader.H>
#include <OrthoFunction.H>
#include <Coefficients.H>
#include <YamlCheck.H>
#include <localmpi.H>
#include <interp.H>

namespace Field
{
  class FieldGenerator;
}

namespace OrthoBasisClasses
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
      Abstract class for any type of basis for evaluating a general
      phase-space function
      
      Provides a factory member that returns derived class instances.
  */
  //! Defines a basis-based 
  class OrthoBasis
  {
    friend class Field::FieldGenerator;
    
  private:

    //! Helper for constructor
    void initialize();

    //! Helper for factory
    // static std::shared_ptr<OrthoBasis>
    // factory_initialize(const YAML::Node& node);

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
    
    //! Total mass on grid
    double totalMass;
    
    //! Evaluate fields in spherical coordinates in centered coordinate system
    virtual std::vector<double>
    sph_eval(double r, double costh, double phi) = 0;


    //! Evaluate fields in cylindrical coordinates in centered coordinate system
    virtual std::vector<double>
    cyl_eval(double r, double costh, double phi) = 0;


    //! Evaluate fields in Cartesian coordinates in centered coordinate system
    virtual std::vector<double>
    crt_eval(double x, double y, double z) = 0;

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

    //! Desired coordinate system for evaluation (default: Spherical)
    enum class Coord { Spherical, Cylindrical, Cartesian, None };

    //! Set the correponding coorindate labels
    static std::map<Coord, std::string> coordLabels;

    //! Get field type
    Coord parseFieldType(std::string Coord_type);
    
    //! Current field coordindate type
    Coord coordinates;

    //! Get field labels
    virtual std::vector<std::string> getFieldLabels(const Coord ctype) = 0;

  public:
    
    //! Constructor from YAML node
    OrthoBasis(const YAML::Node& conf);
    
    //! Constructor from YAML string
    OrthoBasis(const std::string& confstr);
    
    //! Destructor
    virtual ~OrthoBasis(void) {}
    
    //! Set the expansion center
    void setCenter(std::vector<double> center)
    { coefctr = center; }
    
    //! Evaluate basis in desired coordinates
    virtual std::vector<double>
    operator()(double x1, double x2, double x3,
	       const Coord ctype=Coord::Spherical)
    {
      if (ctype==Coord::Spherical)
	return sph_eval(x1, x2, x3);
      else if (ctype==Coord::Cylindrical)
	return cyl_eval(x1, x2, x3);
      else if (ctype==Coord::Cartesian)
	return crt_eval(x1, x2, x3);
      else {
	return sph_eval(x1, x2, x3);
      };
    };
    
    //! Retrieve the coefficients 
    virtual CoefClasses::CoefStrPtr getCoefficients() = 0;

    //! Accumulate new coefficients
    virtual void accumulate(double mass,
			    double x, double y, double z,
			    double u, double v, double w) = 0;
    
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
    static std::shared_ptr<OrthoBasis> factory(const YAML::Node&  conf);
    static std::shared_ptr<OrthoBasis> factory_string(const std::string& conf);
    //@}
    
    //! Set the particle selector callback
    void setSelector(const Callback& f) { ftor = f; }

    //! Clear the particle selector callback
    void clrSelector() { ftor = nullptr; }
  };
  
  /**
     Velocity basis expansion
  */
  class VelocityBasis : public OrthoBasis
  {
    
  public:

    using BasisMap   = std::map<std::string, Eigen::VectorXd>;
    using BasisArray = std::vector<std::vector<BasisMap>>;

  private:
    
    //! Helper for constructor
    void initialize();

    std::shared_ptr<OrthoFunction> ortho;
    std::shared_ptr<Linear1d> interp;
    std::function<double(double)> densfunc;

    std::string model, filename;
    int lmax, nmax, dof;
    double rmin, rmax, ascl, scale, delta;
    
    std::vector<Eigen::VectorXcd> store;
    using dataType = Eigen::Tensor<std::complex<double>, 3>;
    using coefType = Eigen::TensorMap<dataType>;
    std::vector<std::shared_ptr<coefType>> coefs;
    CoefClasses::CoefStrPtr coefstr;

    std::vector<int> usedT;
    std::vector<double> massT;
    int used;
    
  protected:

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "VelocityBasis"; }

    //! Subspace index
    virtual const std::string harmonic()
    { if (dof==3) return "l"; return "m"; }

    virtual std::vector<double>
    sph_eval(double r, double costh, double phi);

    //! Evaluate fields in cylindrical coordinates in centered coordinate system
    virtual std::vector<double>
    cyl_eval(double r, double z, double phi);

    //! Evaluate fields in Cartesian coordinates in centered coordinate system
    virtual std::vector<double>
    crt_eval(double x, double y, double z);

  public:
    
    //! Constructor from YAML node
    VelocityBasis(const YAML::Node& conf);
    
    //! Constructor from YAML string
    VelocityBasis(const std::string& confstr);
    
    //! Destructor
    virtual ~VelocityBasis(void) {}
    
    //! Prescaling factor
    void set_scale(const double scl) { scale = scl; }
    
    //! Zero out coefficients to prepare for a new expansion
    void reset_coefs(void);
    
    //! Make coefficients after accumulation
    void make_coefs(void);
    
    //! Retrieve the coefficients 
    CoefClasses::CoefStrPtr getCoefficients();

    //! Accumulate new coefficients
    virtual void accumulate(double mass,
			    double x, double y, double z,
			    double u, double v, double w);
			    
    //! Get the field labels
    std::vector<std::string> getFieldLabels(const Coord ctype);
  
    //! Return current maximum harmonic order in expansion
    int getLmax() { return lmax; }
    
    //! Return current maximum order in radial expansion
    int getNmax() { return nmax; }

    //! Return potential-density pair of a vector of a vector of 1d
    //! basis-function grids for VelocityBasis, logarithmically spaced
    //! between [logxmin, logxmax] (base 10).
    BasisArray getBasis
    (double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    /*
    std::vector<Eigen::MatrixXd> orthoCheck(int knots=40)
    {
    }
    */

  };
  
  using OrthoBasisPtr = std::shared_ptr<OrthoBasis>;
}
// END: namespace OrthoBasisClasses
  
#endif // _OrthoBasisFactory_H
  
  
  
