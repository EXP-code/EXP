#ifndef _COEFSE_H
#define _COEFSE_H

#include <memory>
#include <Eigen/Eigen>

namespace Coefs
{
  
  //! This is a general struct for any dimensionality
  class CoefStruct
  {
  public:
    //! For YAML config
    std::shared_ptr<char[]> buf;
    
    //! Name of the force
    std::string id;
    
    //! Time stamp
    double time;
    
    typedef
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
    EigenColMajor;
    
    EigenColMajor coefs;
    
    virtual bool read(std::istream& in, bool exp_type) = 0;
  };
  
  //! Specialization of CoefsE for spheres
  class SphStruct : public CoefStruct
  {
  public:
    //! Angular and radial dimension
    int lmax, nmax;
    
    //! Radial scale factor
    double scale;
    
    //! Is basis unit normed?
    bool normed;
    
    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type);
  };
  
  //! Specialization of CoefsE for cylinders
  class CylStruct : public CoefStruct
  {
  public:
    //! Angular and radial dimension
    int mmax, nmax;
    
    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type);
  };
  
  typedef std::shared_ptr<CoefStruct> CoefStrPtr;
  typedef std::shared_ptr<SphStruct>  SphStrPtr;
  typedef std::shared_ptr<CylStruct>  CylStrPtr;
}  
// END namespace Coefs


#endif
  
