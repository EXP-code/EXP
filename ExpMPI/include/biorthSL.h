// This may look like C code, but it is really -*- C++ -*-

#ifndef _biorthSL_h

#define _biorthSL_h 1

#ifdef __GNUG__
#pragma interface
#endif

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

#include <localmpi.h>

#include <numerical.h>
#include <Vector.h>
#include <biorth.h>

#include <sltableMP2.h>


class SphereSL : public AxiSymBiorth
{

private:
  
  int lmax, nmax, numr;
  double rmin, rmax;

  func_1d sphpot, sphdpot, sphdens;

  double xmin, xmax, dxi;

  Vector r;
  Vector xi;
  Vector p0;
  Vector d0;

  struct TableSph* table;

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

  //  void bomb(string oops);

public:

				// Global MPI stuff
  static int mpi;		// default: 0=off

				// Check for cached table
  static int cache;		// default: 1=yes

  static double rscale;		// Radial mapping scale

				// Constructors

  SphereSL(int lmax, int nmax, int numr, double rmin, double rmax,
	   func_1d pot, func_1d dpot, func_1d dens);
  ~SphereSL(void);

				// Required members

  double potl(const int n, const int l, const double x);
  double dens(const int n, const int l, const double x);

  inline double potlR(const int n, const int l, const double r) {
    return potl(n, l, r_to_rb(r) ); }

  inline double potlRZ(const int n, const int l, const double r,
		       const double z) {
    return potl(n, l, r_to_rb(sqrt(r*r+z*z)) ); }


  void potl(const int n, const int l, const double x, Vector& t);
  void dens(const int n, const int l, const double x, Vector& t);

  

  double r_to_rb(double const r);
  double rb_to_r(double const x);
  double d_r_to_rb(double const x);

  double get_potl(const double r, const int l, const Vector& coef);
  double get_dens(const double r, const int l, const Vector& coef);

  inline double rb_min(void) { return xmin; }
  inline double rb_max(void) { return xmax; }

  inline double krnl(const int n, const int l) { return 1.0; }
  inline double norm(const int n, const int l) { return 1.0; }


				// Supplemental functions

  double eigenvalue(int l, int n) {return table[l].ev[n];}
  double d_rb_to_r(const double x);

};

#endif // _BiorthSL_h

