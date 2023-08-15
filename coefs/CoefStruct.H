#ifndef _COEF_STRUCT_H
#define _COEF_STRUCT_H

#include <memory>
#include <complex>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For general rectangular grids

namespace CoefClasses
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
    
    //@{
    //! Coefficient data types
    using EigenColMajor =
      Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    using Eigen3Tensor =
      Eigen::Tensor<std::complex<double>, 3>;
    //@}
    
    //@{
    //! Coefficient data
    EigenColMajor coefs;	// 2d
    Eigen3Tensor  coefT;	// 3d
    //@}
    
    //! Center data
    std::vector<double> ctr;

    //! Destructor
    virtual ~CoefStruct() {}

    //! Read data from a stream
    virtual bool read(std::istream& in, bool exp_type,
		      bool verbose=false) = 0;

    //! Create an empty data storage
    virtual void create() = 0;

    //! Make a complete copy into a new instance
    virtual std::shared_ptr<CoefStruct> deepcopy() = 0;

    //! Copy base-class fields
    void copyfields(std::shared_ptr<CoefStruct> p);

    //! Zero the existing data
    void zerodata() { coefs.setZero(); }

    //! Read-write access to coefficient data (no copy)
    Eigen::Ref<EigenColMajor> setCoefs() { return coefs; }

    //! Read-only access to coefficient data (no copy)
    Eigen::Ref<const EigenColMajor> getCoefs() { return coefs; }

    //! Read-only access to coefficient data (no copy)
    const Eigen3Tensor getCoefT() { return coefT; }

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
    
    //! Constructor
    SphStruct() : lmax(0), nmax(0), scale(1.0), normed(true)
    { geom = "sphere"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

  };
  
  //! Specialization of CoefStruct for cylinders
  class CylStruct : public CoefStruct
  {
  public:
    //! Angular and radial dimension
    int mmax, nmax;
    
    //! Constructor
    CylStruct() : mmax(0), nmax(0)
    { geom = "cylinder"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();
  };
  
  //! Specialization of CoefStruct for slabs
  class SlabStruct : public CoefStruct
  {
  public:
    //! Elements in each dimension
    int nmaxx, nmaxy, nmaxz;
    
    //! Coefficient tensor

    //! Constructor
    SlabStruct() : nmaxx(0), nmaxy(0), nmaxz(0)
    { geom = "slab"; }

    //! No need to read
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();
  };
  
  //! Specialization of CoefStruct for a simple box/cube
  class BoxStruct : public CoefStruct
  {
  public:
    //! Elements in each dimension
    int nmaxx, nmaxy, nmaxz;
    
    //! Coefficient tensor

    //! Constructor
    BoxStruct() : nmaxx(0), nmaxy(0), nmaxz(0)
    { geom = "box"; }

    //! No need to read
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();
  };
  
  //! Specialization of CoefStruct for a data table
  class TblStruct : public CoefStruct
  {
  public:
    //! Number of columns
    int cols;
    
    //! Constructor
    TblStruct() : cols(0) { geom = "table"; }

    //! Read data from a stream
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();
  };
  
  typedef std::shared_ptr<CoefStruct> CoefStrPtr;
  typedef std::shared_ptr<SphStruct>  SphStrPtr;
  typedef std::shared_ptr<CylStruct>  CylStrPtr;
  typedef std::shared_ptr<SlabStruct> SlabStrPtr;
  typedef std::shared_ptr<BoxStruct>  BoxStrPtr;
  typedef std::shared_ptr<TblStruct>  TblStrPtr;

}  
// END namespace CoefClasses


#endif
  
