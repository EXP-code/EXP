// This may look like C code, but it is really -*- C++ -*-

#ifndef _SLGridMP_h

#define _SLGridMP_h 1

#ifdef __GNUG__
#pragma interface
#endif

#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <String.h>
#include <mpi.h>
#include <localmpi.h>

#include <vector.h>

#include <sltableMP.h>

class SLGrid
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

  struct Table** table;

  void init_table(void);
  void compute_table(struct Table* table, int M, int K);
  void compute_table_slave(void);


				// Local MPI stuff
  void mpi_setup(void);
  void mpi_unpack_table(void);
  int mpi_pack_table(struct Table* table, int m, int k);

  int mpi_myid, mpi_numprocs;
  int mpi_bufsz;
  char *mpi_buf;

  void bomb(String oops);

public:

				// Global MPI stuff
  static int mpi;		// default: 0=off

				// Constructors

  SLGrid(int mmax, int nmax, int numr, int numk, double rmin, double rmax,
	 double l);
  ~SLGrid();

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

#endif // _SLGridMP_h
