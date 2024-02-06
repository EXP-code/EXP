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
  
  using VelocityBasisPtr = std::shared_ptr<VelocityBasis>;
}
// END: namespace BasisClasses
  
#endif // _OrthoBasisFactory_H
  
  
  
