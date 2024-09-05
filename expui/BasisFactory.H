#ifndef _BasisFactory_H
#define _BasisFactory_H

#include <functional>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids
#include <yaml-cpp/yaml.h>

#include <ParticleReader.H>
#include <OrthoFunction.H>
#include <Coefficients.H>
#include <YamlCheck.H>
#include <localmpi.H>
#include <exputils.H>
#include <interp.H>

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
      Abstract class for any type of basis for evaluating a general
      phase-space function
      
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
    static std::shared_ptr<Basis>
    factory_initialize(const YAML::Node& node);

  protected:
    
    //! The callback particle selector
    Callback ftor;
    
    //! Return instance for coefficients created from arrays
    CoefClasses::CoefStrPtr coefret;

    //! Current particle index
    unsigned long coefindx;
    
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

    //! Turn on midplane evaluation
    bool midplane = false;

    //! Midplane escursion parameter
    double colh = 4.0;

  public:
    
    //! Constructor from YAML node
    Basis(const YAML::Node& conf,     const std::string& name="Basis");
    
    //! Constructor from YAML string
    Basis(const std::string& confstr, const std::string& name="Basis");
    
    //! Destructor
    virtual ~Basis(void) {}
    
    //! Set the expansion center
    void setCenter(std::vector<double> center)
    { coefctr = center; }
    
    //! Evaluate basis in desired coordinates
    virtual std::vector<double>
    operator()(double x1, double x2, double x3,
	       const Coord ctype=Coord::Spherical);
    
    //! Evaluate fields at a point
    virtual std::vector<double> getFields(double x, double y, double z);
    
    //! Evaluate fields at a point for all coefficients sets
    virtual std::tuple<std::map<std::string, Eigen::VectorXd>,
		       Eigen::VectorXd> getFieldsCoefs
    (double x, double y, double z, std::shared_ptr<CoefClasses::Coefs> coefs);
    
    //! Evaluate fields at a point, and provide field lables
    virtual std::tuple<std::vector<double>, std::vector<std::string>>
    evaluate(double x, double y, double z)
    { return {getFields(x, y, z), getFieldLabels(coordinates)}; }
    
    //! Retrieve the coefficients 
    virtual CoefClasses::CoefStrPtr getCoefficients()
    { return coefret; }

    //! Accumulate new coefficients from full phase space
    virtual void accumulate(double mass,
			    double x, double y, double z,
			    double u, double v, double w) = 0;
    
    //! Accumulate new coefficients from coordinates only
    virtual void accumulate(double x, double y, double z,
			    double mass) = 0;
    
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
    static std::shared_ptr<Basis> factory(const YAML::Node& conf);
    static std::shared_ptr<Basis> factory_string(const std::string& conf);
    //@}
    
    //! Generate coeffients from a particle reader
    virtual CoefClasses::CoefStrPtr
    createFromReader(PR::PRptr reader, std::vector<double> ctr) = 0;

    //! Generate coefficients from a phase-space table
    virtual void
    initFromArray(std::vector<double> ctr) = 0;

    //! Accumulate coefficient contributions from arrays
    virtual void
    addFromArray(Eigen::VectorXd& m, RowMatrixXd& p, bool roundrobin,
		 bool posvelrows) = 0;

    //! Generate coeffients from an array and optional center location
    //! for the expansion
    CoefClasses::CoefStrPtr createFromArray
    (Eigen::VectorXd& m, RowMatrixXd& p, double time=0.0,
     std::vector<double> center={0.0, 0.0, 0.0},
     bool roundrobin=true, bool posvelrows=false);

    //! Create and the coefficients from the array accumulation with the
    //! provided time value
    CoefClasses::CoefStrPtr makeFromArray(double time=0.0);
    
    //! Set the particle selector callback
    void setSelector(const Callback& f) { ftor = f; }

    //! Clear the particle selector callback
    void clrSelector() { ftor = nullptr; }

    //! Turn on/off midplane evaluation (only effective for disk basis)
    void setMidplane(bool value) { midplane = value; }

    //! Height above/below the plane for midplane search in disk scale
    //! lengths
    void setColumnHeight(double value) { colh = value; }
  };
  
  using BasisPtr = std::shared_ptr<Basis>;
}
// END: namespace BasisClasses
  
#endif // _BasisFactory_H
  
  
  
