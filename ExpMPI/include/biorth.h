// This may look like C code, but it is really -*- C++ -*-

// Biorthnormal function class defintion


#ifndef _biorth_h
#define _biorth_h 1

#undef RNUM

static const char rcsid_biorth[] = "$Id$";

#include <string>
#include <Vector.h>
#include <logic.h>


class Biorth
{
public:

  string BiorthID;

  virtual ~Biorth() {}

  virtual double potl(const int, const int, const double) = 0;
  virtual double dens(const int, const int, const double) = 0;
  virtual double krnl(const int, const int) = 0;
  virtual double norm(const int, const int) = 0;
  
				// Report ID

  inline const string& id(void) {return BiorthID;}

				// Error function

  void bomb(const char *s) {
    cerr << "ERROR from " << BiorthID << ": " << s << '\n';
    exit(-1);
  }

};


class AxiSymBiorth : public Biorth
{

friend class BiorthWake;

protected:
  int dof;

public:
  AxiSymBiorth(int DOF=0) { dof = DOF; }

  virtual ~AxiSymBiorth() {}

  virtual double potl(const int, const int, const double) = 0;
  virtual double dens(const int, const int, const double) = 0;
  virtual double potlR(const int , const int, const double) = 0;
  virtual double potlRZ(const int, const int, const double, const double) = 0;
				// Fast for expansion! 
				// (using recursion relation)
  virtual void potl(const int, const int, const double, Vector&) = 0;
  virtual void dens(const int, const int, const double, Vector&) = 0;

  virtual double rb_to_r(double const) = 0;
  virtual double r_to_rb(double const) = 0;
  virtual double d_r_to_rb(double const) = 0;

  virtual double rb_min(void) = 0;
  virtual double rb_max(void) = 0;

  int get_dof(void) { return dof; }

  double get_potl(double r, int l, const Vector& coef);
  double get_dens(double r, int l, const Vector& coef);

};


				// This method tends to be about 3x faster
				// than the direct calls

class BiorthGrid : public AxiSymBiorth
{
private:
  
  AxiSymBiorth* t;

  double rmin;			// Range limits
  double rmax;
  double xmin;
  double xmax;

  int lmax;			// Grid sizes
  int nmax;
  int rnum;

  Matrix* potl_grid;		// Grid storage
  Matrix* potl_grid2;
  Matrix* dens_grid;
  Matrix* dens_grid2;

  Vector x_grid;

  Matrix krnl_grid;
  Matrix norm_grid;		// Stored as the sqr root of normalization
				// for efficiency in wake reconstruction

public:

  BiorthGrid(AxiSymBiorth& T,
	     double RMIN=0.0, double RMAX=20, 
	     int NMAX=10, int LMAX=10, int RNUM=400);
	     
  virtual ~BiorthGrid() {}

				// Required functions

  double potl(const int nn, const int l, const double rb);
  double dens(const int nn, const int l, const double rb);

  inline double potlR(const int nn, const int l, const double r) {
    return potl(nn, l, r_to_rb(r) ); }
  inline double densR(const int nn, const int l, const double r) {
    return dens(nn, l, r_to_rb(r) ); }
  inline double potlRZ(const int nn, const int l, const double r, 
		      const double z) {
    return t->potlRZ(nn, l, r, z); }
  inline void potl(const int nn, const int l, const double r, 
		      Vector& a) {t->potl(nn, l, r, a);}
  inline void dens(const int nn, const int l, const double r, 
		      Vector& a) {t->dens(nn, l, r, a);}

  inline double rb_to_r(double const rb) { return t->rb_to_r(rb); }
  inline double r_to_rb(double const r)  { return t->r_to_rb(r); }
  inline double d_r_to_rb(double const r) { return t->d_r_to_rb(r); }

  inline double rb_min(void) { return xmin; }
  inline double rb_max(void) { return xmax; }

				// Supplemental functions

  inline double krnl(const int nn, const int l) { return krnl_grid[nn][l]; }
  double norm(const int nn, const int l)
    { return norm_grid[nn][l]*norm_grid[nn][l]; }


};




class CBSphere : public AxiSymBiorth
{
private:

  double rbmin;
  double rbmax;

public:

  CBSphere(void);
  virtual ~CBSphere() {}

				// Required functions

  double potl(const int nn, const int l, const double rb);
  double dens(const int nn, const int l, const double rb);
  double krnl(const int nn, const int l);
  double norm(const int nn, const int l);
  void potl(const int nn, const int l, const double r, Vector& a);
  void dens(const int nn, const int l, const double r, Vector& a);

				// Supplemental functions

  inline double potlR(const int nn, const int l, const double r) {
    return potl(nn, l, r_to_rb(r) ); }
  inline double densR(const int nn, const int l, const double r) {
    return dens(nn, l, r_to_rb(r) ); }
  inline double potlRZ(const int nn, const int l, const double r,
		       const double z) {
    return potl(nn, l, r_to_rb(sqrt(r*r+z*z)) ); }

  double rb_to_r(double const);
  double r_to_rb(double const);
  double d_r_to_rb(double const);

  double rb_min(void) {return rbmin;}
  double rb_max(void) {return rbmax;}

};

class HQSphere : public AxiSymBiorth
{
private:

  double rbmin;
  double rbmax;

public:

  HQSphere(void);
  virtual ~HQSphere() {}
				// Required functions

  double potl(const int nn, const int l, const double rb);
  double dens(const int nn, const int l, const double rb);
  double krnl(const int nn, const int l);
  double norm(const int nn, const int l);
  void potl(const int nn, const int l, const double r, Vector& a);
  void dens(const int nn, const int l, const double r, Vector& a);

				// Supplemental functions

  inline double potlR(const int nn, const int l, const double r) {
    return potl(nn, l, r_to_rb(r) ); }
  inline double densR(const int nn, const int l, const double r) {
    return dens(nn, l, r_to_rb(r) ); }
  inline double potlRZ(const int nn, const int l, const double r,
		       const double z) {
    return potl(nn, l, r_to_rb(sqrt(r*r+z*z)) ); }

  double rb_to_r(double const);
  double r_to_rb(double const);
  double d_r_to_rb(double const);

  double rb_min(void) {return rbmin;}
  double rb_max(void) {return rbmax;}

};


class BSSphere : public AxiSymBiorth
{
private:

  double rmax;
  int nmax;
  int lmax;

  Vector *a;
  Matrix t_f, t_y, t_g;
  Vector t_dr;
  int t_n;

  void setup_potl_table(void);

public:

  BSSphere(double RMAX=1.0, int NMAX=10, int LMAX=10);
  virtual ~BSSphere();

				// Required functions

  double potl(const int nn, const int l, const double rb);
  double dens(const int nn, const int l, const double rb);
  inline double krnl(const int nn, const int l) { return 1.0; }
  inline double norm(const int nn, const int l) { return 1.0; }
  void potl(const int nn, const int l, const double r, Vector& t);
  void dens(const int nn, const int l, const double r, Vector& t);

				// Supplemental functions

  inline double potlR(const int nn, const int l, const double r) {
    return potl(nn, l, r ); }
  inline double densR(const int nn, const int l, const double r) {
    return dens(nn, l, r ); }
  inline double potlRZ(const int nn, const int l, const double r,
		       const double z) {
    return potl(nn, l, r_to_rb(sqrt(r*r+z*z)) ); }

  inline double rb_to_r(double const r) { return r; }
  inline double r_to_rb(double const r) { return r; }
  inline double d_r_to_rb(double const r) { return 1.0; }

  inline double rb_min(void) {return 0.0;}
  inline double rb_max(void) {return rmax;}

};

enum ScalarType {density, potential};
string& ScalarTypeName(ScalarType i);

enum BiorthFcts3d {bessel, clutton_brock, hernquist, sturm};

Vector scalar_prod(ScalarType type, double rmin, double rmax, int l, int m,
		   AxiSymBiorth& s, double (*func)(double, int, int), 
		   int nc, int ng);

static string BiorthFcts3dName[] = {"BSSphere", "CBSphere", "HQSphere", "SphereSL"};
		   
#endif
