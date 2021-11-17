// Biorthnormal function class defintion

#ifndef _biorth_H
#define _biorth_H

#undef RNUM

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <memory>
#include <Eigen/Eigen>


class Biorth
{
public:

  std::string BiorthID;

  virtual ~Biorth() {}

  virtual double potl(const int, const int, const double) = 0;
  virtual double dens(const int, const int, const double) = 0;
  virtual double krnl(const int, const int) = 0;
  virtual double norm(const int, const int) = 0;
  
				// Report ID

  inline const std::string& id(void) {return BiorthID;}

				// Error function

  void bomb(const char *s) {
    std::cerr << "ERROR from " << BiorthID << ": " << s << std::endl;
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
  virtual void potl(const int, const int, const double, Eigen::VectorXd&) = 0;
  virtual void dens(const int, const int, const double, Eigen::VectorXd&) = 0;

  virtual double rb_to_r(double const) = 0;
  virtual double r_to_rb(double const) = 0;
  virtual double d_r_to_rb(double const) = 0;

  virtual double rb_min(void) = 0;
  virtual double rb_max(void) = 0;

  int get_dof(void) { return dof; }

  double get_potl(double r, int l, const Eigen::VectorXd& coef);
  double get_dens(double r, int l, const Eigen::VectorXd& coef);

};

typedef std::shared_ptr<AxiSymBiorth> AxiSymBioPtr;

				// This method tends to be about 3x faster
				// than the direct calls

class BiorthGrid : public AxiSymBiorth
{
private:
  
  AxiSymBioPtr t;

  double rmin;			// Range limits
  double rmax;
  double xmin;
  double xmax;

  int lmax;			// Grid sizes
  int nmax;
  int rnum;

  std::vector<Eigen::MatrixXd> potl_grid;		// Grid storage
  std::vector<Eigen::MatrixXd> potl_grid2;
  std::vector<Eigen::MatrixXd> dens_grid;
  std::vector<Eigen::MatrixXd> dens_grid2;

  Eigen::VectorXd x_grid;

  Eigen::MatrixXd krnl_grid;
  Eigen::MatrixXd norm_grid; // Stored as the sqr root of normalization
				// for efficiency in wake reconstruction

public:

  BiorthGrid(AxiSymBioPtr T,
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
		   Eigen::VectorXd& a) {t->potl(nn, l, r, a);}
  inline void dens(const int nn, const int l, const double r, 
		   Eigen::VectorXd& a) {t->dens(nn, l, r, a);}

  inline double rb_to_r(double const rb) { return t->rb_to_r(rb); }
  inline double r_to_rb(double const r)  { return t->r_to_rb(r); }
  inline double d_r_to_rb(double const r) { return t->d_r_to_rb(r); }

  inline double rb_min(void) { return xmin; }
  inline double rb_max(void) { return xmax; }

				// Supplemental functions

  inline double krnl(const int n, const int l) { return krnl_grid(n, l); }
  double norm(const int n, const int l)
  { return norm_grid(n, l)*norm_grid(n, l); }


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
  void potl(const int nn, const int l, const double r, Eigen::VectorXd& a);
  void dens(const int nn, const int l, const double r, Eigen::VectorXd& a);

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
  void potl(const int nn, const int l, const double r, Eigen::VectorXd& a);
  void dens(const int nn, const int l, const double r, Eigen::VectorXd& a);

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

  std::vector<Eigen::VectorXd> a;
  Eigen::MatrixXd t_f, t_y, t_g;
  Eigen::VectorXd t_dr;
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
  void potl(const int nn, const int l, const double r, Eigen::VectorXd& t);
  void dens(const int nn, const int l, const double r, Eigen::VectorXd& t);

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
std::string& ScalarTypeName(ScalarType i);

enum BiorthFcts3d {bessel, clutton_brock, hernquist, sturm};

Eigen::VectorXd scalar_prod(ScalarType type, double rmin, double rmax,
			    int l, int m,
			    AxiSymBiorth& s, double (*func)(double, int, int), 
			    int nc, int ng);

static std::string BiorthFcts3dName[] = {"BSSphere", "CBSphere", "HQSphere", "SphereSL"};
		   
typedef std::shared_ptr<AxiSymBiorth> AxiSymBioPtr;

#endif
