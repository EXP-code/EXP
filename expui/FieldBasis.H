#ifndef _FieldBasis_H
#define _FieldBasis_H

#include <functional>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids
#include <yaml-cpp/yaml.h>

#include "ParticleReader.H"
#include "OrthoFunction.H"
#include "Coefficients.H"
#include "YamlCheck.H"
#include "localmpi.H"
#include "interp.H"

#include "BasisFactory.H"

namespace BasisClasses
{
  /**
     General basis expansion
  */
  class FieldBasis : public Basis
  {
    
  public:

    using BasisArray = Eigen::MatrixXd;
    using PS3 = std::array<double, 3>;
    using PSFunction = std::function<std::vector<double>
				     (double, PS3&, PS3&)>;

  private:
    
    //@{
    //! Helpers for constructor
    void initialize();
    void configure();
    //@}

    std::shared_ptr<OrthoFunction> ortho;
    std::shared_ptr<Linear1d> interp;
    std::function<double(double)> densfunc;

    std::vector<int> usedT;
    std::vector<double> massT;
    int used;
    
  protected:

    PSFunction fieldFunc;
    std::vector<std::string> fieldLabels;
  
    //@{
    //! Parameters
    std::string model, modelname;
    int lmax, nmax;
    double rmin, rmax, ascl, rmapping, delta;
    //@}

    //@{
    //! Coefficient storage
    int nt;
    std::vector<Eigen::VectorXcd> store;
    using dataType = Eigen::Tensor<std::complex<double>, 3>;
    using coefType = Eigen::TensorMap<dataType>;
    std::vector<std::shared_ptr<coefType>> coefs;
    void allocateStore();
    //@}

    //! Number of data fields
    int nfld, dof;

    //! Load coefficients into the new CoefStruct
    virtual void load_coefs(CoefClasses::CoefStrPtr coefs, double time);

    //! Set coefficients
    virtual void set_coefs(CoefClasses::CoefStrPtr coefs);

    //! Valid keys for YAML configurations
    static const std::set<std::string> valid_keys;

    //! Return readable class name
    virtual const std::string classname() { return "FieldBasis"; }

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

    //! Get the field labels
    std::vector<std::string> getFieldLabels(const Coord ctype)
    { return fieldLabels; }
  
  public:
    
   //! Constructor from YAML node
    FieldBasis(const YAML::Node& conf,
	       const std::string& name="FieldBasis") : Basis(conf, name)
    { configure(); }

    //! Constructor from YAML string
    FieldBasis(const std::string& confstr,
	       const std::string& name="FieldBasis") : Basis(confstr, name)
    { configure(); }
    
    //! Destructor
    virtual ~FieldBasis(void) {}
    
    //! Register phase-space functionoid
    void addPSFunction(PSFunction func, std::vector<std::string>& labels);

    //! Coordinate mapping factor
    void set_scale(const double scl) { rmapping = scl; }
    
    //! Zero out coefficients to prepare for a new expansion
    void reset_coefs(void);
    
    //! Make coefficients after accumulation
    void make_coefs(void);
    
    //! Retrieve the coefficients 
    CoefClasses::CoefStrPtr getCoefficients();

    //! Generate coeffients from a particle reader and optional center
    //! location for the expansion
    CoefClasses::CoefStrPtr createFromReader
    (PR::PRptr reader, Eigen::Vector3d center=Eigen::Vector3d::Zero(),
     RowMatrix3d rot=RowMatrix3d::Identity());
    
    
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
    (Eigen::VectorXd& m, RowMatrixXd& p, bool roundrobin=true, bool posvelrows=false);
    
    //! Accumulate new coefficients
    virtual void accumulate(double x, double y, double z,
			    double u, double v, double w,
			    double mass, unsigned long int indx = 0);
			    
    [[deprecated("not relevant for this class")]]
    virtual void accumulate(double x, double y, double z, double mass,
			    unsigned long indx = 0)
    {
    }

    //! Return current maximum harmonic order in expansion
    int getLmax() { return lmax; }
    
    //! Return current maximum order in radial expansion
    int getNmax() { return nmax; }

    //! Return potential-density pair of a vector of a vector of 1d
    //! basis-function grids for FieldBasis, logarithmically spaced
    //! between [logxmin, logxmax] (base 10).
    BasisArray getBasis
    (double logxmin=-3.0, double logxmax=0.5, int numgrid=2000);

    //! Compute the orthogonality of the basis by returning inner
    //! produce matrices
    Eigen::MatrixXd orthoCheck() { return ortho->testOrtho(); }
  };
  
  class VelocityBasis : public FieldBasis
  {

  private:

    void assignFunc();

  public:
    //! Constructor from YAML node
    VelocityBasis(const YAML::Node& conf);

    //! Constructor from YAML string
    VelocityBasis(const std::string& confstr);
  };

  using FieldBasisPtr    = std::shared_ptr<FieldBasis>;
  using VelocityBasisPtr = std::shared_ptr<VelocityBasis>;
}
// END: namespace BasisClasses
  
#endif // _FieldBasis_H
  
  
  
