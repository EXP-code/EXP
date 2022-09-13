#ifndef _COEF_STRUCT_H
#define _COEF_STRUCT_H

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
    
    //! Geometry
    std::string geom;
    
    //! Name of the force
    std::string id;
    
    //! Time stamp
    double time;
    
    typedef
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
    EigenColMajor;
    
    //! Coefficient data
    EigenColMajor coefs;
    
    //! Center data
    std::vector<double> ctr;

    //! Read data from a stream
    virtual bool read(std::istream& in, bool exp_type) = 0;

    //! Copy base-class fields
    void copyfields(std::shared_ptr<CoefStruct> p);

  };
  
  //! Specialization of CoefStruct for spheres
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

    //! Copy all of data contents to a new instance
    std::shared_ptr<SphStruct> deepcopy();

  };
  
  //! Specialization of CoefStruct for cylinders
  class CylStruct : public CoefStruct
  {
  public:
    //! Angular and radial dimension
    int mmax, nmax;
    
    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type);

    //! Copy all of data contents to a new instance
    std::shared_ptr<CylStruct> deepcopy();
  };
  
  //! Specialization of CoefStruct for a data table
  class TblStruct : public CoefStruct
  {
  public:
    //! Number of columns
    int cols;
    
    //! No need to read
    bool read(std::istream& in, bool exp_type);

    //! Copy all of data contents to a new instance
    std::shared_ptr<TblStruct> deepcopy();
  };
  
  typedef std::shared_ptr<CoefStruct> CoefStrPtr;
  typedef std::shared_ptr<SphStruct>  SphStrPtr;
  typedef std::shared_ptr<CylStruct>  CylStrPtr;
  typedef std::shared_ptr<TblStruct>  TblStrPtr;

}  
// END namespace Coefs


#endif
  
