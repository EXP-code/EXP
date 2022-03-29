#ifndef _Bessel_H
#define _Bessel_H

#include <assert.h>
#include <SphericalBasis.H>

double *sbessjz(int n, int m);
double sbessj(int n, double x);

class MixtureBasis;

/** Computes the potential, acceleration and density using eigenfunctions of the Spherical Laplacian: spherical Bessel functions
*/
class Bessel : public SphericalBasis
{

private:

  void get_pot_coefs(int l, double *coef, double *p, double *dp);
  void get_pot_coefs_safe(int l, double *coef, double *p, double *dp,
			      double **potd1, double **dpot1);
  void get_dens_coefs(int l, double *coef, double *p);
  void get_dpotl(int lmax, int nmax, double r, Eigen::MatrixXd& p, Eigen::MatrixXd& dp, int tid);
  void get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);
  void get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid);
  void get_potl_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, Eigen::MatrixXd& d,int tid);
  double get_dens(double r, int l, double *coef);

  void initialize();
  //  double norm(int, int);
  //  double knl(int, int);

  bool firstime_coef;
  bool firstime_accel;

  //! Grid to hold tabulated basis
  class RGrid {
  public:
    Eigen::MatrixXd rw;
    Eigen::MatrixXd rw2;
    int nmax;
  };
  
  std::vector<RGrid> dens_grid, potl_grid;
  Eigen::VectorXd r_grid;
  double r_grid_del;

  //! Cache roots for spherical Bessel functions
  class Roots {
  public:

    int l;
    int n;
    double *a;

    Roots(int L, int nmax) : l(L), n(nmax) {
      a = sbessjz(l-1,n);
    }

    ~Roots() { delete [] (a+1);}
  };

  Roots *p;

  double dens(double r, int n);
  double potl(double r, int n);

  void make_grid(double rmin, double rmax, int lmax, int nmax);
  
public:
  
  //! Number of entries in fixed table (static variable)
  static int RNUM;

  //! Constructor
  Bessel(Component* c0, const YAML::Node& conf, MixtureBasis* m=0);

  //! Destructor
  virtual ~Bessel();
};

#endif // Bessel.H
