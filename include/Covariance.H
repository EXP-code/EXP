#ifndef _Covariance_H
#define _Covariance_H

#include <functional>
#include <tuple>
#include <map>
#include <set>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids
#include <yaml-cpp/yaml.h>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>
#include <highfive/eigen.hpp>

#include "localmpi.H"
#include "exputils.H"

namespace BasisClasses
{
  /** 
      Class to maintain subsample covariance of basis functions
  */
  class SubsampleCovariance
  {
  public:
    
    using ParamCallback = std::function<void (HighFive::File& file)>;
    
    using CoefCovarType = std::tuple<Eigen::VectorXcd, Eigen::MatrixXcd>;
    
    using CovarData = std::tuple<Eigen::VectorXi, Eigen::VectorXd,
				 std::vector<std::vector<CoefCovarType>>>;
    
  protected:
    
    ParamCallback paramCallback;
    
    //! Variance file versioning (Version 1.0 has more granularity and
    //! poor compression)
    inline static const std::string CovarianceFileVersion = "1.1";
    
    //! Store covariance matrices?
    bool covar = true;
    
    //! Store summed covariance only?
    bool summed = true;

    //! Use 32-bit float storage for covariance?
    bool floatType = false;
    
    //! Fixed-point multiplier for time mapping.  Eight decimal places
    //! should be enough here...
    const double multiplier = 1.0e+08;
    
    //! Round time key to emulated fixed-point arithmetic
    double roundTime(double time)
    {
      return std::floor(time * multiplier + 0.5) / multiplier;
    }
    
    //! Round time key to emulated fixed-point arithmetic
    int roundTimeInt(double time)
    {
      return std::floor(time * multiplier + 0.5);
    }
    
    //@{
    //! HDF5 compression settings
    unsigned H5compress  = 5;
    unsigned H5chunk     = 1024*1024; // (1 MB)
    bool     H5shuffle   = true;
    bool     H5szip      = false;
    //@}
    
    //! Basis identifier string
    std::string BasisID;
    
    //! Time list
    std::vector<double> times;
    
    //! Imported data members
    std::map<double, CovarData> covarData;
    
    //! Write H5 covariance; returns updated group count
    virtual unsigned writeCovarH5
    (HighFive::Group& group, CovarData& elem, unsigned count, double time);
    
    //! Write H5 parameters for coefficient covariance
    virtual void writeCovarH5Params(HighFive::File& file)
    {
      if (paramCallback) {
	paramCallback(file);
      } else {
	throw std::runtime_error("SubsampleCovariance::writeCovarH5Params: no paramCallback defined");
      }
    }
    
    //! Add coefficient covariance data to an HDF5 file
    virtual void extendCoefCovariance
    (const std::string& filename, CovarData& elem, double time);
    
  public:
    
    //! Constructor from scratch
    SubsampleCovariance(ParamCallback func, const std::string& BasisID,
			bool flt, bool sum=false, bool cov=false);
    
    //! Constructor from HDF5 file
    SubsampleCovariance(const std::string& filename, int stride=1);
    
    //! Get time list
    std::vector<double> Times() { return times; }
    
    //! Return a vector of tuples of basis functions and the
    //! covariance matrix for subsamples of particles
    virtual std::vector<std::vector<CoefCovarType>>
    getCoefCovariance(double time)
    {
      auto it = covarData.find(roundTime(time));
      if (it == covarData.end())
	throw std::runtime_error("SubsampleCovariance::getCoefCovariance: "
				 "time not found");
      return std::get<2>(it->second);
    }
    
    //! Get sample counts for the covariance computation
    virtual std::tuple<Eigen::VectorXi, Eigen::VectorXd>
    getCovarSamples(double time)
    {
      auto it = covarData.find(roundTime(time));
      if (it == covarData.end())
	throw std::runtime_error("SubsampleCovariance::getCovarSamples: "
				 "time not found");
      return std::make_tuple(std::get<0>(it->second), std::get<1>(it->second));
    }
    
    //! Write coefficient covariance data to an HDF5 file, creating a
    //! filename based on component name and run tag
    void writeCoefCovariance(const std::string& compname,
			     const std::string& runtag,
			     CovarData& covarElem,
			     double time=0.0);
    
    //! Write coefficient covariance data to an HDF5 file with a
    //! supplied path + filename
    void writeCoefCovariance(const std::string& filename,
			     CovarData& covarElem,
			     double time=0.0);

    //! Choose between float and double storage for covariance
    void setCovarFloatType(bool flt) { floatType = flt; }
    
    //! HDF5 compression settings
    void setCovarH5Compress(unsigned level, unsigned chunksize, bool shuffle, bool szip=false)
    {
      H5compress = level;
      H5chunk    = chunksize;
      H5shuffle  = shuffle;
      H5szip     = szip;
    }

    //! Return the basis ID string that generated this covariance sample
    const std::string& basisIDname() const { return BasisID; }
  };
}
// END: namespace BasisClasses
  
#endif // _Covariance_H
  
