// This may look like C code, but it is really -*- C++ -*-

#ifndef _SLGridMP_h

#define _SLGridMP_h 1

#ifdef __GNUG__
#pragma interface
#endif

#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <string>
#include <mpi.h>
#include <localmpi.h>

#include <Vector.h>
// #include <Matrix.h>

#include <massmodel.h>
#include <sltableMP2.h>

class SLGridCyl
{

private:
  
  int mmax, nmax, numr, numk;
  double rmin, rmax, l;

  double dk;

  double xmin, xmax, dxi;

  Vector kv;
  Vector r;
  Vector xi;
  Vector p0;
  Vector d0;

  struct TableCyl** table;

  void init_table(void);
  void compute_table(struct TableCyl* table, int M, int K);
  void compute_table_slave(void);
  int read_cached_table(void);
  void write_cached_table(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int mpi_pack_table(struct TableCyl* table, int m, int k);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  char *mpi_buf;

  void bomb(string oops);

public:

				// Global MPI stuff
  static int mpi;		// default: 0=off

				// Check for cached table
  static int cache;		// default: 1=yes

				// Exponential scale length
  static double A;		// default: 1.0

				// Coordinate map
  static int cmap;		// default: 1=on

				// Constructors

  SLGridCyl(int mmax, int nmax, int numr, int numk, double rmin, double rmax,
	 double l);
  ~SLGridCyl();

				// Members

  double eigenvalue(int m, int k, int n) {return table[m][k].ev[n];}
  double r_to_xi(double r);
  double xi_to_r(double x);
  double d_xi_to_r(double x);

  double get_pot(double x, int m, int n, int k, int which=1);
  double get_dens(double x, int m, int n, int k, int which=1);
  double get_force(double x, int m, int n, int k, int which=1);

  void get_pot(Vector& vec, double x, int m, int k, int which=1);
  void get_dens(Vector& vec, double x, int m, int k, int which=1);
  void get_force(Vector& vec, double x, int m, int k, int which=1);

  void get_pot(Matrix& tab, double x, int m, int which=1);
  void get_dens(Matrix& tab, double x, int m, int which=1);
  void get_force(Matrix& tab, double x, int m, int which=1);

  void get_pot(Matrix* tab, double x, int mMin, int mMax, int which=1);
  void get_dens(Matrix* tab, double x, int mMin, int mMax, int which=1);
  void get_force(Matrix* tab, double x, int mMin, int mMax, int which=1);

};


class SLGridSph
{

private:
  
  int lmax, nmax, numr;
  double rmin, rmax;

  double xmin, xmax, dxi;

  Vector r;
  Vector xi;
  Vector p0;
  Vector d0;

  struct TableSph* table;

  void initialize(int LMAX, int NMAX, int NUMR,
		  double RMIN, double RMAX);

  void init_table(void);
  void compute_table(struct TableSph* table, int L);
  void compute_table_slave(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int mpi_pack_table(struct TableSph* table, int l);
  int read_cached_table(void);
  void write_cached_table(void);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  char *mpi_buf;

  void bomb(string oops);

public:

				// Global MPI stuff
  static int mpi;		// default: 0=off

				// Check for cached table
  static int cache;		// default: 1=yes

				// Coordinate map
  static int cmap;		// default: 0=off

				// Constructors

  SLGridSph(int lmax, int nmax, int numr, double rmin, double rmax,
	    SphericalModelTable *mod);
  SLGridSph(int lmax, int nmax, int numr, double rmin, double rmax);
  ~SLGridSph();

				// Members

  double eigenvalue(int l, int n) {return table[l].ev[n];}
  double r_to_xi(double r);
  double xi_to_r(double x);
  double d_xi_to_r(double x);

  double get_pot(double x, int l, int n, int which=1);
  double get_dens(double x, int l, int n, int which=1);
  double get_force(double x, int l, int n, int which=1);

  void get_pot(Vector& vec, double x, int l, int which=1);
  void get_dens(Vector& vec, double x, int l, int which=1);
  void get_force(Vector& vec, double x, int l, int which=1);

  void get_pot(Matrix& tab, double x, int which=1);
  void get_dens(Matrix& tab, double x, int which=1);
  void get_force(Matrix& tab, double x, int which=1);

};


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

  struct TableSlab** table;

  void init_table(void);
  void compute_table(struct TableSlab* table, int kx, int ky);
  void compute_table_slave(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int mpi_pack_table(struct TableSlab* table, int kx, int ky);
  int read_cached_table(void);
  void write_cached_table(void);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  char *mpi_buf;

  void bomb(string oops);

public:

				// Global MPI stuff
  static int mpi;		// default: 0=off

				// Check for cached table
  static int cache;		// default: 1=yes

  static double H;		// Scale height, default=0.1

  static double L;		// Periodic box size, default=1.0

  static double ZBEG;		// Offset from origin, default=1.0e-4

  static double ZEND;		// Potential offset, default=1.0e-5


				// Constructors

  SLGridSlab(int kmax, int nmax, int numz, double zmax);
  ~SLGridSlab();

				// Members

  double eigenvalue(int kx, int ky, int n) {return table[kx][ky].ev[n];}
  double z_to_xi(double z);
  double xi_to_z(double x);
  double d_xi_to_z(double x);

  double get_pot(double x, int kx, int ky, int n, int which=1);
  double get_dens(double x, int kx, int ky, int n, int which=1);
  double get_force(double x, int kx, int ky, int n, int which=1);

  void get_pot(Vector& vec, double x, int kx, int ky, int which=1);
  void get_dens(Vector& vec, double x, int kx, int ky, int which=1);
  void get_force(Vector& vec, double x, int kx, int ky, int which=1);

  void get_pot(Matrix& tab, double x, int which=1);
  void get_dens(Matrix& tab, double x, int which=1);
  void get_force(Matrix& tab, double x, int which=1);

};


#endif // _SLGridMP_h

