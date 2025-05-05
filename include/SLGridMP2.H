#ifndef _SLGridMP_h
#define _SLGridMP_h

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

#include <mpi.h>
#include <localmpi.H>

#include <config_exp.h>

#include <massmodel.H>
#include <sltableMP2.H>
#include <yaml-cpp/yaml.h>

#include <libvars.H>
using namespace __EXP__;

#if HAVE_LIBCUDA==1
#include <cudaUtil.cuH>
#include <cudaMappingConstants.cuH>
#endif


//!! Spherical SL grid class
class SLGridSph
{

private:
  
  int lmax, nmax, numr;
  double rmin, rmax;

  int cmap, diverge;
  double rmap, dfac;

  double xmin, xmax, dxi;

  Eigen::VectorXd r;
  Eigen::VectorXd xi;
  Eigen::VectorXd p0;
  Eigen::VectorXd d0;

  // Use shared pointers to avoid copying and memory leaks
  using table_ptr_1D = std::shared_ptr<TableSph[]>;
  table_ptr_1D table;

  void initialize(int LMAX, int NMAX, int NUMR,
		  double RMIN, double RMAX, 
		  bool CACHE, int CMAP, double RMAP);

  void init_table(void);
  void compute_table(TableSph* table, int L);
  void compute_table_worker(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int  mpi_pack_table(TableSph* table, int l);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  std::shared_ptr<char []> mpi_buf;

  void bomb(string oops);

  //! Use basis cache
  bool cache;

  //! For deep debugging
  bool tbdbg;

  //! Default model file name
  const std::string default_model = "SLGridSph.model";

  //! Default cache file name
  const std::string default_cache = ".slgrid_sph_cache";

  //! Model file name
  std::string model_file_name;

  //! Cache file name
  std::string sph_cache_name;

  //! Write HDF5 cache
  void WriteH5Cache();

  //! Read HDF5 cache
  bool ReadH5Cache();

  //! Cache versioning
  inline static const std::string Version = "1.0";

public:

  //! Flag for MPI enabled (default: 0=off)
  static int mpi;

				// Constructors

  //! Constructor with model table
  SLGridSph(std::shared_ptr<SphericalModelTable> mod,
	    int lmax, int nmax, int numr, double rmin, double rmax,
	    bool cache, int Cmap, double RMAP,
	    std::string cachename=".slgrid_sph_cache",
	    bool Verbose=false);

  //! Constructor (uses file *model_file_name* for file)
  SLGridSph(std::string modelname,
	    int lmax, int nmax, int numr, double rmin, double rmax,
	    bool cache, int Cmap, double RMAP, 
	    int DIVERGE, double DFAC,
	    std::string cachename=".slgrid_sph_cache",
	    bool Verbose=false);

  //! Constructor (from cache file)
  SLGridSph(std::string cachename);


  //! Destructor
  virtual ~SLGridSph();

  //! Read the cache and report parameters
  static YAML::Node getHeader(const std::string& cachefile);

				// Members

  //! Eigenvalue for index and harmonic order l
  double eigenvalue(int l, int n);

  //! Map radial coordinate to nondimensional coordinate
  double r_to_xi(double r);

  //! Map nondimensional coordinate to radial coordinate
  double xi_to_r(double x);

  //! Jacobian of nondimensional coordinate mapping
  double d_xi_to_r(double x);

  
  //! Get potential for dimensionless coord with harmonic order l and radial orer n
  double get_pot(double x, int l, int n, int which=1);

  //! Get density for dimensionless coord with harmonic order l and radial orer n  
  double get_dens(double x, int l, int n, int which=1);

  //! Get force for dimensionless coord with harmonic order l and radial orer n
  double get_force(double x, int l, int n, int which=1);

  /** Get potential for dimensionless coord with harmonic order l and radial orer n
      Return all radial order values in Vector
  */
  void get_pot(Eigen::VectorXd& vec, double x, int l, int which=1);
  /** Get density for dimensionless coord with harmonic order l and radial orer n
      Return all radial order values in Vector
  */

  //! Read and print the cache and return the header parameters as a
  //! map/dictionary
  static std::map<std::string, std::string>
  cacheInfo(const std::string& cachefile, bool verbose=true);

  //! Compute the orthogonality of the basis by returning inner
  //! produce matrices
  std::vector<Eigen::MatrixXd> orthoCheck(int knots=40);

#if HAVE_LIBCUDA==1
  void initialize_cuda(std::vector<cudaArray_t>& cuArray,
		       thrust::host_vector<cudaTextureObject_t>& tex);

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    cudaMappingConstants ret;

    ret.rscale = rmap;
    ret.hscale = 0.0;
    ret.xmin   = xmin;
    ret.xmax   = xmax;
    ret.ymin   = 0.0;
    ret.ymax   = 0.0;
    ret.numr   = numr;
    ret.numx   = 0;
    ret.numy   = 0;
    ret.dxi    = dxi;
    ret.dyi    = 0.0;
    ret.cmapR  = cmap;

    return ret;
  }

#endif

  void get_dens(Eigen::VectorXd& vec, double x, int l, int which=1);
  /** Get force for dimensionless coord with harmonic order l and radial orer n
      Return all radial order values in Vector
  */
  void get_force(Eigen::VectorXd& vec, double x, int l, int which=1);

  /** Get potential for dimensionless coord with harmonic order l and radial order n
      Return Matrix with first dim harmonic order and second dim radial order
  */
  void get_pot(Eigen::MatrixXd& tab, double x, int which=1);

  /** Get density for dimensionless coord with harmonic order l and radial order n
      Return Matrix with first dim harmonic order and second dim radial order
  */
  void get_dens(Eigen::MatrixXd& tab, double x, int which=1);

  /** Get force for dimensionless coord with harmonic order l and radial order n
      Return Matrix with first dim harmonic order and second dim radial order
  */
  void get_force(Eigen::MatrixXd& tab, double x, int which=1);

  //@{
  //! Get the current minimum and maximum radii for the expansion
  double getRmin() { return rmin; }
  double getRmax() { return rmax; }
  //@}

  //@{
  //! Get expansion limits
  int getLmax() { return lmax; }
  int getNmax() { return nmax; }
  //@}
};


//! Target density models for slabs
class SlabModel
{
protected:
  //! For cache identification
  std::string id;

public:

  //! Factory constructor
  static std::shared_ptr<SlabModel> createModel(const std::string type);

  //! Return the potential
  virtual double pot(double z) = 0;

  //! Return the derivative of the potential
  virtual double dpot(double z) = 0;

  //! Return the surface density
  virtual double dens(double z) = 0;

  //! Get model ID
  std::string ID() const { return id; }
};

//! Slab (one-dimensional) SL grid class
class SLGridSlab
{

private:
  
  int numk, nmax, numz;
  double zmax;

  double xmin, xmax, dxi;

  Eigen::VectorXd z;
  Eigen::VectorXd xi;
  Eigen::VectorXd p0;
  Eigen::VectorXd d0;

  // Use shared pointers to avoid copying and memory leaks
  using table_ptr_1D = std::shared_ptr<TableSlab[]>;
  using table_ptr_2D = std::shared_ptr<table_ptr_1D[]>;

  table_ptr_2D table;

  void init_table(void);
  void compute_table(TableSlab* table, int kx, int ky);
  void compute_table_worker(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int  mpi_pack_table(TableSlab* table, int kx, int ky);
  bool ReadH5Cache(void);
  void WriteH5Cache(void);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;

  std::shared_ptr<char[]> mpi_buf;

  void bomb(string oops);

  //! For deep debugging
  bool tbdbg;

  //@{

  //! This defines a class providing selectable maps from the infinite
  //! to finite interval
  enum CoordMapTypes {Tanh, Sech, Linear};

  //! The base class for all coordinate maps
  class CoordMap
  {
  protected:
    //! Scale parameter
    double H;

  public:
    //! Constructor
    CoordMap(double H) : H(H) {}
  
    //! Destructor
    virtual ~CoordMap() {}

    //! Convert from vertical to mapped coordinate
    virtual double z_to_xi  (double z) = 0;

    //! Convert from mapped coordinate back to vertical
    virtual double xi_to_z  (double z) = 0;

    //! Jacobian of the transformation
    virtual double d_xi_to_z(double z) = 0;

    //! Coordinate map factory
    static std::unique_ptr<CoordMap> factory
    (CoordMapTypes type, double H);
  };


  //! x = tanh(z/H)
  class TanhMap : public CoordMap
  {
  public:
    TanhMap(double H) : CoordMap(H) {}
    virtual double z_to_xi  (double z);
    virtual double xi_to_z  (double z);
    virtual double d_xi_to_z(double z);
  };

  //! x = z/sqrt(z^2 + H^2)
  class SechMap : public CoordMap
  {
  public:
    SechMap(double H) : CoordMap(H) {}
    virtual double z_to_xi  (double z);
    virtual double xi_to_z  (double z);
    virtual double d_xi_to_z(double z);
  };

  //! x = z
  class LinearMap : public CoordMap
  {
  public:
    LinearMap(double H) : CoordMap(H) {}
    virtual double z_to_xi  (double z);
    virtual double xi_to_z  (double z);
    virtual double d_xi_to_z(double z);
  };

  //! Model type
  std::string type;

  //! The current map created by the SLGridSlab constructor
  std::unique_ptr<CoordMap> mM;

public:

  //! Global MPI flag, default: 0=off
  static int mpi;	

  //! Check for cached table, default: 1=yes
  static int cache;		

  //! Scale height, default=0.1
  static double H;

  //! Periodic box size, default=1.0
  static double L;

  //! Offset from origin, default=1.0e-4
  static double ZBEG;

  //! Potential offset, default=1.0e-5
  static double ZEND;


  //! Constructor
  SLGridSlab(int kmax, int nmax, int numz, double zmax,
	     const std::string type="isothermal", bool Verbose=false);

  //! Destructor
  ~SLGridSlab();


  //@{
  //! Members

  //! Get n^th eigenvalue for given wave number indices
  double eigenvalue(int kx, int ky, int n);

  //! Map from vertical coordinate to dimensionless coordinate
  double z_to_xi(double z)   { return mM->z_to_xi(z); }

  //! Map from dimensionless coordinate to vertical coordinate
  double xi_to_z(double x)   { return mM->xi_to_z(x); }

  //! Jacobian of coordinate mapping
  double d_xi_to_z(double x) { return mM->d_xi_to_z(x); }

  //! Get potential for dimensionless coord with given wave numbers and index
  double get_pot(double x, int kx, int ky, int n, int which=1);

  //! Get density for dimensionless coord with given wave numbers and index
  double get_dens(double x, int kx, int ky, int n, int which=1);

  //! Get force for dimensionless coord with given wave numbers and index
  double get_force(double x, int kx, int ky, int n, int which=1);

  /** Get potential for member for dimensionless coord with given wave numbers
      Return Vector for all indices
  */
  void get_pot(Eigen::VectorXd& vec, double x, int kx, int ky, int which=1);
  /** Get density for dimensionless coord with given wave numbers
      Return Vector for all indices
  */
  void get_dens(Eigen::VectorXd& vec, double x, int kx, int ky, int which=1);
  /** Get force for dimensionless coord with given wave numbers
      Return Vector for all indices
  */
  void get_force(Eigen::VectorXd& vec, double x, int kx, int ky, int which=1);

  /** Get potential for dimensionless coord with given wave numbers
      Return Matrix with first dimension containing x and y wavenumbers 
      packed with y index varying most quicly, second index is vertical order.
  */
  void get_pot(Eigen::MatrixXd& tab, double x, int which=1);
  /** Get density for dimensionless coord with given wave numbers
      Return Matrix with first dimension containing x and y wavenumbers 
      packed with y index varying most quicly, second index is vertical order.
  */
  void get_dens(Eigen::MatrixXd& tab, double x, int which=1);
  /** Get force for dimensionless coord with given wave numbers
      Return Matrix with first dimension containing x and y wavenumbers 
      packed with y index varying most quicly, second index is vertical order.
  */
  void get_force(Eigen::MatrixXd& tab, double x, int which=1);

  //! Compute the orthogonality of the basis by returning inner
  //! produce matrices
  std::vector<Eigen::MatrixXd> orthoCheck(int knots=40);

  //@}

#if HAVE_LIBCUDA==1
  void initialize_cuda(std::vector<cudaArray_t>& cuArray,
		       thrust::host_vector<cudaTextureObject_t>& tex);

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    cudaMappingConstants ret;

    ret.hscale = H;
    ret.xmin   = xmin;
    ret.xmax   = xmax;
    ret.ymin   = 0.0;
    ret.ymax   = 0.0;
    ret.numr   = numz;
    ret.numx   = 0;
    ret.numy   = 0;
    ret.dxi    = dxi;
    ret.dyi    = 0.0;

    return ret;
  }

#endif

};


#endif // _SLGridMP_h

