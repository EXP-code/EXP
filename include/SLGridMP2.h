// This may look like C code, but it is really -*- C++ -*-

#ifndef _SLGridMP_h
#define _SLGridMP_h 1

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

#include <mpi.h>
#include <localmpi.h>

#include <Vector.h>
// #include <Matrix.h>

#include <massmodel.h>
#include <sltableMP2.h>

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

  Vector kv;
  Vector r;
  Vector xi;
  Vector p0;
  Vector d0;

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

public:

  //! Global MPI indicator, default: 0=off
  static int mpi;

  //! Check for cached table, default: 1=yes
  static int cache;		

  //! Exponential scale length, default: 1.0
  static double A;		


  //! Constructor
  SLGridCyl(int mmax, int nmax, int numr, int numk, double rmin, double rmax,
	 double l, int Cmap=0, double Scale=1.0);
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
  void get_pot(Vector& vec, double x, int m, int k, int which=1);
  //! Fill vector with desired density basis 
  void get_dens(Vector& vec, double x, int m, int k, int which=1);
  //! Fill vector with desired force basis 
  void get_force(Vector& vec, double x, int m, int k, int which=1);

  //! Fill Matrix with desired potential basis 
  void get_pot(Matrix& tab, double x, int m, int which=1);
  //! Fill Matrix with desired potential basis 
  void get_dens(Matrix& tab, double x, int m, int which=1);
  //! Fill Matrix with desired potential basis 
  void get_force(Matrix& tab, double x, int m, int which=1);

  //! Fill Matricies with desired potential basis 
  void get_pot(Matrix* tab, double x, int mMin, int mMax, int which=1);
  //! Fill Matricies with desired density basis 
  void get_dens(Matrix* tab, double x, int mMin, int mMax, int which=1);
  //! Fill Matricies with desired force basis 
  void get_force(Matrix* tab, double x, int mMin, int mMax, int which=1);

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

  Vector r;
  Vector xi;
  Vector p0;
  Vector d0;

  TableSph* table;

  void initialize(int LMAX, int NMAX, int NUMR,
		  double RMIN, double RMAX, int CMAP, double SCALE);

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

  bool my_model;

  void bomb(string oops);

public:

  //! Flag for MPI enabled (default: 0=off)
  static int mpi;

  //! Check for cached table (default: 1=yes)
  static int cache;

  //! Model file name
  static string model_file_name;

  //! Cache file name
  static string sph_cache_name;

				// Constructors

  //! Constructor with model table
  SLGridSph(int lmax, int nmax, int numr, double rmin, double rmax,
	    SphericalModelTable *mod, int Cmap=0, double Scale=1.0);
  //! Constructor (uses file *model_file_name* for file)
  SLGridSph(int lmax, int nmax, int numr, double rmin, double rmax,
	    int Cmap=0, double Scale=1.0, int DIVERGE=0, double DFAC=1.0);
  //! Destructor
  ~SLGridSph();

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
  void get_pot(Vector& vec, double x, int l, int which=1);
  /** Get density for dimensionless coord with harmonic order l and radial orer n
      Return all radial order values in Vector
  */
  void get_dens(Vector& vec, double x, int l, int which=1);
  /** Get force for dimensionless coord with harmonic order l and radial orer n
      Return all radial order values in Vector
  */
  void get_force(Vector& vec, double x, int l, int which=1);

  /** Get potential for dimensionless coord with harmonic order l and radial order n
      Return Matrix with first dim harmonic order and second dim radial order
  */
  void get_pot(Matrix& tab, double x, int which=1);

  /** Get density for dimensionless coord with harmonic order l and radial order n
      Return Matrix with first dim harmonic order and second dim radial order
  */
  void get_dens(Matrix& tab, double x, int which=1);

  /** Get force for dimensionless coord with harmonic order l and radial order n
      Return Matrix with first dim harmonic order and second dim radial order
  */
  void get_force(Matrix& tab, double x, int which=1);

};


//! Slab (one-dimensional) SL grid class
class SLGridSlab
{

private:
  
  int numk, nmax, numz;
  double zmax;

  double xmin, xmax, dxi;

  Vector z;
  Vector xi;
  Vector p0;
  Vector d0;

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
  SLGridSlab(int kmax, int nmax, int numz, double zmax);

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
  void get_pot(Vector& vec, double x, int kx, int ky, int which=1);
  /** Get density for dimensionless coord with given wave numbers
      Return Vector for all indices
  */
  void get_dens(Vector& vec, double x, int kx, int ky, int which=1);
  /** Get force for dimensionless coord with given wave numbers
      Return Vector for all indices
  */
  void get_force(Vector& vec, double x, int kx, int ky, int which=1);

  /** Get potential for dimensionless coord with given wave numbers
      Return Matrix with first dimension containing x and y wavenumbers 
      packed with y index varying most quicly, second index is vertical order.
  */
  void get_pot(Matrix& tab, double x, int which=1);
  /** Get density for dimensionless coord with given wave numbers
      Return Matrix with first dimension containing x and y wavenumbers 
      packed with y index varying most quicly, second index is vertical order.
  */
  void get_dens(Matrix& tab, double x, int which=1);
  /** Get force for dimensionless coord with given wave numbers
      Return Matrix with first dimension containing x and y wavenumbers 
      packed with y index varying most quicly, second index is vertical order.
  */
  void get_force(Matrix& tab, double x, int which=1);

};


#endif // _SLGridMP_h

