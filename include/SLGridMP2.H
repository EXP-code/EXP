#ifndef _SLGridMP_h
#define _SLGridMP_h

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

#include <mpi.h>
#include <localmpi.H>

#include <config.h>

#include <massmodel.H>
#include <sltableMP2.H>

#if HAVE_LIBCUDA==1
#include <cudaUtil.cuH>
#include <cudaMappingConstants.cuH>
#endif

//! Cylindrical SL grid class
class SLGridCyl
{

private:
  
  int mmax, nmax, numr, numk;
  double rmin, rmax, l;

  int cmap;
  double scale;

  double dk;

  double xmin, xmax, dxi;

  Eigen::VectorXd kv;
  Eigen::VectorXd r;
  Eigen::VectorXd xi;
  Eigen::VectorXd p0;
  Eigen::VectorXd d0;

  TableCyl** table;

  void init_table(void);
  void compute_table(TableCyl* table, int M, int K);
  void compute_table_slave(void);
  int read_cached_table(void);
  void write_cached_table(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int mpi_pack_table(TableCyl* table, int m, int k);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  char *mpi_buf;

  void bomb(string oops);

  bool cache;		
  bool tbdbg;

public:

  //! Global MPI indicator, default: 0=off
  static int mpi;

  //! Exponential scale length, default: 1.0
  static double A;		


  //! Constructor
  SLGridCyl(int mmax, int nmax, int numr, int numk, double rmin, double rmax,
	    double l, bool cache, int Cmap, double Scale,
	    bool Verbose=false);
  //! Destructor
  ~SLGridCyl();

				// Members

  //! Return eigenvalue for given order
  double eigenvalue(int m, int k, int n) {return table[m][k].ev[n];}
  //! Dimensional to dimensionless coordinate mapping
  double r_to_xi(double r);
  //! Dimensionless to dimensional coordinate mapping
  double xi_to_r(double x);
  //! Jacobian of dimensionless mapping
  double d_xi_to_r(double x);

  //! Get potential basis value
  double get_pot(double x, int m, int n, int k, int which=1);
  //! Get density basis value
  double get_dens(double x, int m, int n, int k, int which=1);
  //! Get force basis value
  double get_force(double x, int m, int n, int k, int which=1);

  //! Fill vector with desired potential basis 
  void get_pot(Eigen::VectorXd& vec, double x, int m, int k, int which=1);
  //! Fill vector with desired density basis 
  void get_dens(Eigen::VectorXd& vec, double x, int m, int k, int which=1);
  //! Fill vector with desired force basis 
  void get_force(Eigen::VectorXd& vec, double x, int m, int k, int which=1);

  //! Fill Matrix with desired potential basis 
  void get_pot(Eigen::MatrixXd& tab, double x, int m, int which=1);

#if HAVE_LIBCUDA==1
  void get_pot_cuda(float* tab, double x, int m, int which=1);
#endif

  //! Fill Matrix with desired potential basis 
  void get_dens(Eigen::MatrixXd& tab, double x, int m, int which=1);
  //! Fill Matrix with desired potential basis 
  void get_force(Eigen::MatrixXd& tab, double x, int m, int which=1);

  //! Fill Matricies with desired potential basis 
  void get_pot(Eigen::MatrixXd* tab, double x, int mMin, int mMax, int which=1);
  //! Fill Matricies with desired density basis 
  void get_dens(Eigen::MatrixXd* tab, double x, int mMin, int mMax, int which=1);
  //! Fill Matricies with desired force basis 
  void get_force(Eigen::MatrixXd* tab, double x, int mMin, int mMax, int which=1);

};


//!! Spherical SL grid class
class SLGridSph
{

private:
  
  int lmax, nmax, numr;
  double rmin, rmax;

  int cmap;
  double scale;

  double xmin, xmax, dxi;

  Eigen::VectorXd r;
  Eigen::VectorXd xi;
  Eigen::VectorXd p0;
  Eigen::VectorXd d0;

  TableSph* table;

  void initialize(int LMAX, int NMAX, int NUMR,
		  double RMIN, double RMAX, 
		  bool CACHE, int CMAP, double SCALE);

  void init_table(void);
  void compute_table(TableSph* table, int L);
  void compute_table_slave(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int mpi_pack_table(TableSph* table, int l);
  int read_cached_table(void);
  void write_cached_table(void);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  char *mpi_buf;

  void bomb(string oops);

  bool cache;
  bool tbdbg;

public:

  //! Flag for MPI enabled (default: 0=off)
  static int mpi;

  //! Model file name
  static string model_file_name;

  //! Cache file name
  static string sph_cache_name;

				// Constructors

  //! Constructor with model table
  SLGridSph(std::shared_ptr<SphericalModelTable> mod,
	    int lmax, int nmax, int numr, double rmin, double rmax,
	    bool cache, int Cmap, double Scale, bool Verbose=false);

  //! Constructor (uses file *model_file_name* for file)
  SLGridSph(int lmax, int nmax, int numr, double rmin, double rmax,
	    bool cache, int Cmap, double Scale, 
	    int DIVERGE, double DFAC, bool Verbose=false);

  //! Destructor
  virtual ~SLGridSph();

				// Members

  //! Eigenvalue for index and harmonic order l
  double eigenvalue(int l, int n) {return table[l].ev[n];}
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

#if HAVE_LIBCUDA==1
  void initialize_cuda(std::vector<cudaArray_t>& cuArray,
		       thrust::host_vector<cudaTextureObject_t>& tex);

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    cudaMappingConstants ret;

    ret.rscale = scale;
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

  TableSlab** table;

  void init_table(void);
  void compute_table(TableSlab* table, int kx, int ky);
  void compute_table_slave(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int mpi_pack_table(TableSlab* table, int kx, int ky);
  int read_cached_table(void);
  void write_cached_table(void);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  char *mpi_buf;

  void bomb(string oops);

  bool tbdbg;

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
  SLGridSlab(int kmax, int nmax, int numz, double zmax, bool Verbose=false);

  //! Destructor
  ~SLGridSlab();

				// Members

  //! Get n^th eigenvalue for given wave number indices
  double eigenvalue(int kx, int ky, int n) {return table[kx][ky].ev[n];}

  //! Map from vertical coordinate to dimensionless coordinate
  double z_to_xi(double z);

  //! Map from dimensionless coordinate to vertical coordinate
  double xi_to_z(double x);

  //! Jacobian of coordinate mapping
  double d_xi_to_z(double x);

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

};


#endif // _SLGridMP_h

