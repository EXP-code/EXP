#ifndef _COEF_STRUCT_H
#define _COEF_STRUCT_H

#include <memory>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids

namespace CoefClasses
{
  //! This is a general struct for any dimensionality
  class CoefStruct
  {
  public:
    //! For YAML config
    std::string buf;
    
    //! Geometry
    std::string geom;
    
    //! Name of the force
    std::string id;
    
    //! Time stamp
    double time;
    
    //! Coefficient data
    Eigen::VectorXcd store;
    
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
    virtual void zerodata() { store.setZero(); }

    //! Read-write access to coefficient data (no copy)
    Eigen::Ref<Eigen::VectorXcd> setCoefs() { return store; }

    //! Set new coefficient data (no copy)
    void setCoefs(Eigen::VectorXcd& STORE)
    {
      if (STORE.size() != store.size())
	throw std::invalid_argument("CoefStruct::setCoefs: coefficient vector size does not match");
      store = STORE;
    }

    //! Read-only access to coefficient data (no copy)
    Eigen::Ref<const Eigen::VectorXcd> getCoefs() { return store; }

  };
  
  //! Specialization of CoefStruct for spheres
  class SphStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::MatrixXcd>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

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

    //! Allocate storage arrays
    void allocate()
    {
      int rows = (lmax+1)*(lmax+2)/2;
      int cols = nmax;
      store.resize(rows*cols);
      coefs = std::make_shared<coefType>(store.data(), rows, cols);
    }

    //! Assign matrix
    void assign(const Eigen::MatrixXcd& mat, int Lmax, int Nmax)
    {
      lmax = Lmax;
      nmax = Nmax;
      allocate();
      *coefs = mat;
    }

  };
  
  //! Specialization of CoefStruct for cylinders
  class CylStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::MatrixXcd>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

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

    //! Allocate storage arrays
    void allocate()
    {
      int rows = mmax+1;
      int cols = nmax;
      store.resize(rows*cols);
      coefs = std::make_shared<coefType>(store.data(), rows, cols);
    }

    //! Assign matrix
    void assign(const Eigen::MatrixXcd& mat, int Mmax, int Nmax)
    {
      mmax = Mmax;
      nmax = Nmax;
      allocate();
      *coefs = mat;
    }
  };
  
  //! Specialization of CoefStruct for a slab
  class SlabStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 3>>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

    //! Angular and radial dimension
    int nmaxx, nmaxy, nmaxz;
    
    //! Basis dimensions
    int nx, ny, nz, dim;

    //! Constructor
    SlabStruct() : nmaxx(0), nmaxy(0), nmaxz(0)
    { geom = "slab"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Allocate storage arrays
    void allocate()
    {
      nx = 2*nmaxx + 1;
      ny = 2*nmaxy + 1;
      nz = nmaxz;
      dim = nx * ny * nz;
      store.resize(dim);
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
    }


    //! Assign matrix
    void assign(const Eigen::Tensor<std::complex<double>, 3>& dat)
    {
      const auto& d = dat.dimensions();
      nx = d[0];
      ny = d[1];
      nz = d[2];

      nmaxx = (nx - 1)/2;
      nmaxy = (ny - 1)/2;
      nmaxz = nz;

      allocate();
      *coefs = dat;
    }

  };
  
  //! Specialization of CoefStruct for a cube
  class CubeStruct : public CoefStruct
  {
  private:
    //! Array dimensions
    int nx, ny, nz, dim;

  public:
    //! Coefficient map type
    using coefType = Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 3>>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

    //! Angular and radial dimension
    int nmaxx, nmaxy, nmaxz;
    
    //! Constructor
    CubeStruct() : nmaxx(0), nmaxy(0), nmaxz(0), nx(0), ny(0), nz(0), dim(0)
    { geom = "slab"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Allocate storage arrays
    void allocate()
    {
      // Grid dimenions
      nx = 2*nmaxx + 1;
      ny = 2*nmaxy + 1;
      nz = 2*nmaxz + 1;
      dim = nx * ny * nz;

      // Allocate
      store.resize(dim);
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
    }

    //! Assign matrix
    void assign(const Eigen::Tensor<std::complex<double>, 3>& dat)
    {
      const auto& d = dat.dimensions();

      // Assign grid dimensions
      nx = d[0];
      ny = d[1];
      nz = d[2];
      dim = nx * ny * nz;

      // Deduce max wave numbers from grid dimensions
      nmaxx = (d[0]-1)/2;
      nmaxy = (d[1]-1)/2;
      nmaxz = (d[2]-1)/2;

      // Allocate storage and assign
      allocate();
      *coefs = dat;
    }
  };
  

  //! Specialization of CoefStruct for a data table
  class TblStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::VectorXcd>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

    //! Number of columns
    int cols;
    
    //! Constructor
    TblStruct() : cols(0) { geom = "table"; }

    //! No need to read
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Assign array
    void assign(const Eigen::VectorXcd& arr)
    {
      cols  = arr.size();
      store = arr;
      coefs = std::make_shared<coefType>(store.data(), store.size());
    }

    //! Allocate storage arrays
    void allocate()
    {
      store.resize(cols);
      coefs = std::make_shared<coefType>(store.data(), cols);
    }

  };
  
  //! Specialization of CoefStruct for collection of trajectories
  class TrajStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::MatrixXcd>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

    //! Number of trajectories
    int traj;
    
    //! Phase-space rank
    int rank;
    
    //! Constructor
    TrajStruct() : traj(0), rank(0) { geom = "trajectory"; }

    //! No need to read
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Assign array
    void assign(const Eigen::MatrixXd& arr)
    {
      traj = arr.rows();
      rank = arr.cols();

      allocate();
      *coefs = arr;
    }

    //! Allocate storage arrays
    void allocate()
    {
      store.resize(traj*rank);
      coefs = std::make_shared<coefType>(store.data(), traj, rank);
    }

  };
  
  //! Specialization of CoefStruct for spheres for general fields
  class SphFldStruct : public CoefStruct
  {
  public:
    //! Coefficient data type
    using dataType = Eigen::Tensor<std::complex<double>, 3>;

    //! Coefficient map type
    using coefType = Eigen::TensorMap<dataType>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

    //! Number of data fields
    int nfld;

    //! Angular and radial dimension
    int lmax, nmax;
    
    //! Radial scale factor
    double scale;
    
    //! Constructor
    SphFldStruct() : nfld(0), lmax(0), nmax(0), scale(1.0)
    { geom = "sphere"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Allocate storage arrays
    void allocate()
    {
      int rows = (lmax+1)*(lmax+2)/2;
      int cols = nmax;
      store.resize(nfld*rows*cols);
      coefs = std::make_shared<coefType>(store.data(), nfld, rows, cols);
    }

    //! Assign data
    void assign(const Eigen::VectorXcd& data, int Nfld, int Lmax, int Nmax)
    {
      nfld = Nfld;
      lmax = Lmax;
      nmax = Nmax;
      allocate();
      store = data;
    }

  };
  

  //! Specialization of CoefStruct for disks for general fields
  class CylFldStruct : public CoefStruct
  {
  public:
    //! Coefficient data type
    using dataType = Eigen::Tensor<std::complex<double>, 3>;

    //! Coefficient map type
    using coefType = Eigen::TensorMap<dataType>;

    //! Coefficient map
    std::shared_ptr<coefType> coefs;

    //! Number of data fields
    int nfld;

    //! Angular and radial dimension
    int mmax, nmax;
    
    //! Radial scale factor
    double scale;
    
    //! Constructor
    CylFldStruct() : nfld(0), mmax(0), nmax(0), scale(1.0)
    { geom = "cylinder"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Allocate storage arrays
    void allocate()
    {
      int rows = mmax+1;
      int cols = nmax;
      store.resize(nfld*rows*cols);
      coefs = std::make_shared<coefType>(store.data(), nfld, rows, cols);
    }

    //! Assign matrix
    void assign(const Eigen::VectorXcd& data, int Nfld, int Mmax, int Nmax)
    {
      nfld = Nfld;
      mmax = Mmax;
      nmax = Nmax;
      allocate();
      store = data;
    }

  };
  

  using CoefStrPtr = std::shared_ptr<CoefStruct>   ;
  using SphStrPtr  = std::shared_ptr<SphStruct>    ;       
  using CylStrPtr  = std::shared_ptr<CylStruct>    ;   
  using SlabStrPtr = std::shared_ptr<SlabStruct>   ;  
  using CubeStrPtr = std::shared_ptr<CubeStruct>   ;  
  using TblStrPtr  = std::shared_ptr<TblStruct>    ;  
  using TrajStrPtr = std::shared_ptr<TrajStruct>   ;  
  using SphFldPtr  = std::shared_ptr<SphFldStruct> ;  
  using PlrFldPtr  = std::shared_ptr<CylFldStruct> ;  

}  
// END namespace CoefClasses


#endif
  
