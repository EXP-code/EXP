#ifndef _COEF_STRUCT_H
#define _COEF_STRUCT_H

#include <memory>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids

namespace CoefClasses
{
  using RowMatrix3d = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;

  class Coefs;

  //! This is a general struct for any dimensionality
  class CoefStruct
  {
  protected:

    //! Coefs instance
    Coefs* C = nullptr;

    //! Default gravitational contant 
    double G = 1.0;

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
    Eigen::Vector3d ctr = Eigen::Vector3d::Zero();

    //! Orientation data
    RowMatrix3d rot = RowMatrix3d::Identity();

    //! Constructor
    CoefStruct(Coefs* C=nullptr) : C(C) {}

    //! Destructor
    virtual ~CoefStruct() {}

    //! Set owner
    void setOwner(Coefs* Cptr) { C = Cptr; }

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

    //! Set new center data
    void setCenter(Eigen::Vector3d& CTR) {ctr = CTR; }

    //! Read-only access to center
    Eigen::Vector3d getCenter() { return ctr; }

    //! Set new rotation matrix
    void setRotation(RowMatrix3d& ROT) { rot = ROT; }

    //! Read-only access to orientation
    RowMatrix3d getRotation() { return rot; }

    //! Set coefficient time (no copy)
    void setTime(double& STORE)
    {
      time = STORE;
    }

    //! Read-only access to center (no copy)
    double getTime() { return time; }

    //! Set default gravitational constant
    void setGravConstant(double val) { G = val; }

    //! Get gravitational constant
    double getGravConstant();

  };
  
  //! Specialization of CoefStruct for spheres
  class SphStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::MatrixXcd>;

    //! Angular and radial dimension
    int lmax, nmax;
    
    //! Radial scale factor
    double scale;
    
    //! Is basis unit normed?
    bool normed;

  private:
    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<coefType> coefs;

  public:
    //! Constructor
    SphStruct(Coefs* C=nullptr) :
      lmax(0), nmax(0), scale(1.0), normed(true), CoefStruct(C)
    { geom = "sphere"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Lazy-create and cache the coefficient map
    const coefType& getCoefMap() const
    {
      if (!coefs) {
        int rows = (lmax+1)*(lmax+2)/2;
        int cols = nmax;
        coefs = std::make_shared<coefType>
          (const_cast<std::complex<double>*>(store.data()), rows, cols);
      }
      return *coefs;
    }

    //! Non-const version for assignment
    coefType& getCoefMapMutable()
    {
      if (!coefs) {
        int rows = (lmax+1)*(lmax+2)/2;
        int cols = nmax;
        coefs = std::make_shared<coefType>(store.data(), rows, cols);
      }
      return *coefs;
    }

    //! Allocate storage arrays
    void allocate()
    {
      int rows = (lmax+1)*(lmax+2)/2;
      int cols = nmax;
      coefs.reset(); // Clear map BEFORE resize
      store.resize(rows*cols);
    }

    //! Assign matrix
    void assign(const Eigen::MatrixXcd& mat, int Lmax, int Nmax)
    {
      lmax = Lmax;
      nmax = Nmax;
      allocate();
      getCoefMapMutable() = mat;
    }

    //! Initialize coefficient map with dimensions
    void initCoefMap(int rows, int cols)
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), rows, cols);
    }

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::MatrixXcd& mat)
    {
      int rows = (lmax+1)*(lmax+2)/2;
      int cols = nmax;
      if (mat.rows() != rows || mat.cols() != cols)
	throw std::invalid_argument("setCoefs: size mismatch");

      getCoefMapMutable() = mat;
    }

    //@{
    //! Direct access to underlying Eigen Map (for client code)
    const coefType& getCoefs() const { return getCoefMap(); }
    coefType& getCoefs() { return getCoefMapMutable(); }
    //@}

  };

  //! Specialization of CoefStruct for cylinders
  class CylStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::MatrixXcd>;

    //! Angular and radial dimension
    int mmax, nmax;

  private:
    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<coefType> coefs;

  public:
    //! Constructor
    CylStruct(Coefs* C=nullptr) : mmax(0), nmax(0), CoefStruct(C)
    { geom = "cylinder"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Lazy-create and cache the coefficient map
    const coefType& getCoefMap() const
    {
      if (!coefs) {
	int rows = mmax+1;
	int cols = nmax;
	coefs->resize(rows, cols);
      }
      return *coefs;
    }

    //! Non-const version for assignment
    coefType& getCoefMapMutable()
    {
      int rows = mmax+1;
      int cols = nmax;
      coefs.reset(); // Clear old map before creating new one
      coefs = std::make_shared<coefType>(store.data(), rows, cols);
      return *coefs;
    }

    //! Allocate storage arrays
    void allocate()
    {
      int rows = mmax+1;
      int cols = nmax;
      
      coefs.reset(); // Clear map BEFORE resize
      store.resize(rows*cols);
    }

    //! Assign matrix
    void assign(const Eigen::MatrixXcd& mat, int Mmax, int Nmax)
    {
      mmax = Mmax;
      nmax = Nmax;
      allocate();
      getCoefMapMutable() = mat;
    }

    //! Initialize coefficient map with dimensions
    void initCoefMap(int rows, int cols)
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), rows, cols);
    }

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::MatrixXcd& mat)
    {
      int rows = mmax+1;
      int cols = nmax;

      if (mat.rows() != rows || mat.cols() != cols)
	throw std::invalid_argument("setCoefs: size mismatch");
      getCoefMapMutable() = mat;
    }

    //@{
    //! Direct access to underlying Eigen Map (for client code)
    const coefType& getCoefs() const { return getCoefMap(); }
    coefType& getCoefs() { return getCoefMapMutable(); }
    //@}

  };
  
  //! Specialization of CoefStruct for a slab
  class SlabStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 3>>;

    //! Angular and radial dimension
    int nmaxx, nmaxy, nmaxz;
    
    //! Basis dimensions
    int nx, ny, nz, dim;

  private:
    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<coefType> coefs;

  public:
    //! Constructor
    SlabStruct(Coefs* C=nullptr) : nmaxx(0), nmaxy(0), nmaxz(0), CoefStruct(C)
    { geom = "slab"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();


// For Map-based classes:
/**/


    //! Lazy-create and cache the coefficient tensor map
    const coefType& getCoefTensorMap() const
    {
      if (!coefs) {
	coefs = std::make_shared<coefType>
	  (const_cast<std::complex<double>*>(store.data()), nx, ny, nz);
      }
      return *coefs;
    }


    //! Non-const version for assignment
    coefType& getCoefTensorMapMutable()
    {
      coefs.reset(); // Clear old map before creating new one
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
      return *coefs;
    }

    //! Allocate storage arrays
    void allocate()
    {
      nx = 2*nmaxx + 1;
      ny = 2*nmaxy + 1;
      nz = nmaxz;
      dim = nx * ny * nz;
      coefs.reset(); // Clear map BEFORE resize
      store.resize(dim);
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
      getCoefTensorMapMutable() = dat;
    }

    //! Initialize coefficient tensor map with dimensions
    void initCoefMap(int nx, int ny, int nz)
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
    }

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::Tensor<std::complex<double>, 3>& tensor)
    {
      const auto& d = tensor.dimensions();
      if (d[0] != nx || d[1] != ny || d[2] != nz)
	throw std::invalid_argument("setCoefs: dimension mismatch");
      getCoefTensorMapMutable() = tensor;
    }

    //! Direct access to underlying Eigen TensorMap (for client code)
    const coefType& getCoefs() const { return getCoefTensorMap(); }
    coefType& getCoefs() { return getCoefTensorMapMutable(); }
  };
  
  //! Specialization of CoefStruct for a cube
  class CubeStruct : public CoefStruct
  {
  private:
    //! Array dimensions
    int nx, ny, nz, dim;

    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 3>>> coefs;

  public:
    //! Coefficient map type
    using coefType = Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 3>>;

    //! Angular and radial dimension
    int nmaxx, nmaxy, nmaxz;

    //! Constructor
    CubeStruct(Coefs* C=nullptr) :
      nmaxx(0), nmaxy(0), nmaxz(0), nx(0), ny(0), nz(0), dim(0), CoefStruct(C)
    { geom = "slab"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Lazy-create and cache the coefficient tensor map
    const coefType& getCoefTensorMap() const
    {
      if (!coefs) {
	coefs = std::make_shared<coefType>
	  (const_cast<std::complex<double>*>(store.data()), nx, ny, nz);
      }
      return *coefs;
    }

    //! Non-const version for assignment
    coefType& getCoefTensorMapMutable()
    {
      coefs.reset(); // Clear old map before creating new one
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
      return *coefs;
    }
    
    //! Allocate storage arrays
    void allocate()
    {
      // Grid dimenions
      nx = 2*nmaxx + 1;
      ny = 2*nmaxy + 1;
      nz = 2*nmaxz + 1;
      dim = nx * ny * nz;

      // Allocate
      coefs.reset(); // Clear map BEFORE resize
      store.resize(dim);
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
      getCoefTensorMapMutable() = dat;
    }

    //! Initialize coefficient tensor map with dimensions
    void initCoefMap(int nx, int ny, int nz)
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), nx, ny, nz);
    }

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::Tensor<std::complex<double>, 3>& tensor)
    {
      const auto& d = tensor.dimensions();
      if (d[0] != nx || d[1] != ny || d[2] != nz)
	throw std::invalid_argument("setCoefs: dimension mismatch");
      getCoefTensorMapMutable() = tensor;
    }

    //! Direct access to underlying Eigen TensorMap (for client code)
    const coefType& getCoefs() const { return getCoefTensorMap(); }
    coefType& getCoefs() { return getCoefTensorMapMutable(); }
  };

  //! Specialization of CoefStruct for a data table
  class TblStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::VectorXcd>;

    //! Number of columns
    int cols;

  private:
    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<coefType> coefs;

  public:
    //! Constructor
    TblStruct(Coefs* C=nullptr) : cols(0), CoefStruct(C) { geom = "table"; }

    //! No need to read
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Lazy-create and cache the coefficient map
    const coefType& getCoefMap() const
    {
      if (!coefs) {
	coefs = std::make_shared<coefType>
	  (const_cast<std::complex<double>*>(store.data()), store.size());
      }
      return *coefs;
    }

    //! Non-const version for assignment
    coefType& getCoefMapMutable()
    {
      coefs.reset(); // Clear old map before creating new one
      coefs = std::make_shared<coefType>(store.data(), store.size());
      return *coefs;
    }

    //! Assign array
    void assign(const Eigen::VectorXcd& arr)
    {
      cols = arr.size();
      coefs.reset(); // Clear map BEFORE resize
      store = arr;
    }

    //! Allocate storage arrays
    void allocate()
    {
      coefs.reset(); // Clear map BEFORE resize
      store.resize(cols);
    }

    //! Initialize coefficient map
    void initCoefMap()
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), store.size());
    }

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::VectorXcd& vec)
    {
      if (vec.size() != cols)
	throw std::invalid_argument("setCoefs: size mismatch");
      getCoefMapMutable() = vec;
    }

    //@{
    //! Direct access to underlying Eigen Map (for client code)
    const coefType& getCoefs() const { return getCoefMap(); }
    coefType& getCoefs() { return getCoefMapMutable(); }
    //@}

  };

  //! Specialization of CoefStruct for collection of trajectories
  class TrajStruct : public CoefStruct
  {
  public:
    //! Coefficient map type
    using coefType = Eigen::Map<Eigen::MatrixXcd>;

    //! Number of trajectories
    int traj;
    
    //! Phase-space rank
    int rank;

private:
    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<coefType> coefs;

  public:
    //! Constructor
    TrajStruct(Coefs* C=nullptr) : traj(0), rank(0), CoefStruct(C)
    { geom = "trajectory"; }

    //! No need to read
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Lazy-create and cache the coefficient map
    const coefType& getCoefMap() const
    {
      if (!coefs) {
	coefs = std::make_shared<coefType>
	  (const_cast<std::complex<double>*>(store.data()), traj, rank);
      }
      return *coefs;
    }

    //! Non-const version for assignment
    coefType& getCoefMapMutable()
    {
      coefs.reset(); // Clear old map before creating new one
      coefs = std::make_shared<coefType>(store.data(), traj, rank);
      
      return *coefs;
    }

    //! Assign array
    void assign(const Eigen::MatrixXd& arr)
    {
      traj = arr.rows();
      rank = arr.cols();
      
      allocate();
      getCoefMapMutable() = arr;
    }

    //! Allocate storage arrays
    void allocate()
    {
      coefs.reset(); // Clear map BEFORE resize
      store.resize(traj*rank);
    }

    //! Initialize coefficient map with dimensions
    void initCoefMap(int rows, int cols)
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), rows, cols);
    }

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::MatrixXcd& mat)
    {
      traj = mat.rows();
      rank = mat.cols();
      
      allocate();
      getCoefMapMutable() = mat;
    }

    //@{
    //! Direct access to underlying Eigen Map (for client code)
    const coefType& getCoefs() const { return getCoefMap(); }
    coefType& getCoefs() { return getCoefMapMutable(); }
    //@}

  };
  
  //! Specialization of CoefStruct for spheres for general fields
  class SphFldStruct : public CoefStruct
  {
  public:
    //! Coefficient data type
    using dataType = Eigen::Tensor<std::complex<double>, 3>;

    //! Coefficient map type
    using coefType = Eigen::TensorMap<dataType>;

    //! Number of data fields
    int nfld;
    //! Angular and radial dimension
    int lmax, nmax;
    
    //! Radial scale factor
    double scale;

  private:
    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<coefType> coefs;

  public:
    //! Constructor
    SphFldStruct(Coefs* C=nullptr) :
      nfld(0), lmax(0), nmax(0), scale(1.0), CoefStruct(C)
    { geom = "sphere"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Lazy-create and cache the coefficient tensor map
    const coefType& getCoefTensorMap() const
    {
      if (!coefs) {
	int rows = (lmax+1)*(lmax+2)/2;
	int cols = nmax;
	coefs = std::make_shared<coefType>
	  (const_cast<std::complex<double>*>(store.data()), nfld, rows, cols);
      }
      return *coefs;
    }

    //! Non-const version for assignment
    coefType& getCoefTensorMapMutable()
    {
      int rows = (lmax+1)*(lmax+2)/2;
      int cols = nmax;
      coefs.reset(); // Clear old map before creating new one
      coefs = std::make_shared<coefType>(store.data(), nfld, rows, cols);
      return *coefs;
    }

    //! Allocate storage arrays
    void allocate()
    {
      int rows = (lmax+1)*(lmax+2)/2;
      int cols = nmax;
      coefs.reset(); // Clear map BEFORE resize
      store.resize(nfld*rows*cols);
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

    //! Initialize coefficient tensor map with dimensions
    void initCoefMap(int nfld, int rows, int cols)
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), nfld, rows, cols);
    }

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::Tensor<std::complex<double>, 3>& tensor)
    {
      int rows = (lmax+1)*(lmax+2)/2;
      int cols = nmax;

      const auto& d = tensor.dimensions();
      if (d[0] != nfld || d[1] != rows || d[2] != cols)
	throw std::invalid_argument("setCoefs: dimension mismatch");

      getCoefTensorMapMutable() = tensor;
    }

    //! Direct access to underlying Eigen TensorMap (for client code)
    const coefType& getCoefs() const { return getCoefTensorMap(); }
    coefType& getCoefs() { return getCoefTensorMapMutable(); }
  };

  //! Specialization of CoefStruct for disks for general fields
  class CylFldStruct : public CoefStruct
  {
  public:
    //! Coefficient data type
    using dataType = Eigen::Tensor<std::complex<double>, 3>;

    //! Coefficient map type
    using coefType = Eigen::TensorMap<dataType>;

    //! Number of data fields
    int nfld;

    //! Angular and radial dimension
    int mmax, nmax;
    
    //! Radial scale factor
    double scale;

  private:
    //! Cached coefficient map (mutable for lazy initialization)
    mutable std::shared_ptr<coefType> coefs;

  public:
    //! Constructor
    CylFldStruct(Coefs* C=nullptr) :
      nfld(0), mmax(0), nmax(0), scale(1.0), CoefStruct(C)
    { geom = "cylinder"; }

    //! Read native EXP coefficients
    bool read(std::istream& in, bool exp_type, bool verbose=false);

    //! Create an empty data storage
    void create();

    //! Copy all of data contents to a new instance
    std::shared_ptr<CoefStruct> deepcopy();

    //! Lazy-create and cache the coefficient tensor map
    const coefType& getCoefTensorMap() const
    {
      if (!coefs) {
	coefs = std::make_shared<coefType>
	  (const_cast<std::complex<double>*>(store.data()), nfld, mmax+1, nmax);
      }
      return *coefs;
    }

    //! Non-const version for assignment
    coefType& getCoefTensorMapMutable()
    {
      int rows = mmax+1;
      int cols = nmax;
      coefs.reset(); // Clear old map before creating new one
      coefs = std::make_shared<coefType>(store.data(), nfld, rows, cols);
      return *coefs;
    }

    //! Allocate storage arrays
    void allocate()
    {
      int rows = mmax+1;
      int cols = nmax;
      coefs.reset(); // Clear map BEFORE resize
      store.resize(nfld*rows*cols);
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

    //! Initialize coefficient tensor map with dimensions
    void initCoefMap(int nfld, int rows, int cols)
    {
      coefs.reset();
      coefs = std::make_shared<coefType>(store.data(), nfld, rows, cols);
    }    

    //! Set the underlying coefficient data
    void setCoefs(const Eigen::Tensor<std::complex<double>, 3>& tensor)
    {
      int rows = mmax+1;
      int cols = nmax;

      const auto& d = tensor.dimensions();
      if (d[0] != nfld || d[1] != rows || d[2] != cols)
	throw std::invalid_argument("setCoefs: dimension mismatch");

      getCoefTensorMapMutable() = tensor;
    }

    //! Direct access to underlying Eigen TensorMap (for client code)
    const coefType& getCoefs() const { return getCoefTensorMap(); }
    coefType& getCoefs() { return getCoefTensorMapMutable(); }
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
