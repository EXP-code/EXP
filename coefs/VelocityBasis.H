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

#include <BasisFactory.H>

namespace BasisClasses
{
  /**
     Velocity basis expansion
  */
  class VelocityBasis : public Basis
  {
    
  public:

    using BasisArray = Eigen::MatrixXd;

  private:
    
    //@{
    //! Helpers for constructor
    void initialize();
    void configure();
    //@}

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
    VelocityBasis(const YAML::Node& conf) : Basis(conf) { configure(); }

    VelocityBasis(const std::string& confstr) : Basis(confstr) { configure(); }
    
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

    //! Generate coeffients from a particle reader and optional center
    //! location for the expansion
    CoefClasses::CoefStrPtr createFromReader
    (PR::PRptr reader, std::vector<double> center={0.0, 0.0, 0.0});
    
    
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
    
    //! Accumulate new coefficients
    virtual void accumulate(double mass,
			    double x, double y, double z,
			    double u, double v, double w);
			    
    [[deprecated("not relevant for this class")]]
    virtual void accumulate(double x, double y, double z, double mass)
    {
    }

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
    Eigen::MatrixXd orthoCheck() { return ortho->testOrtho(); }
  };
  
  using VelocityBasisPtr = std::shared_ptr<VelocityBasis>;
}
// END: namespace BasisClasses
  
#endif // _VelocityBasisFactory_H
  
  
  
