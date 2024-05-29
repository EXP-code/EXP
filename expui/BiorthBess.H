#ifndef _BiorthBess_H_
#define _BiorthBess_H_

#include <vector>
#include <memory>
#include <cmath>

#include <Eigen/Eigen>

Eigen::VectorXd sbessjz(int n, int m);

class BiorthBess
{
private:

  //! Grid to hold tabulated basis
  class RGrid
  {
  public:
    Eigen::MatrixXd rw;
    Eigen::MatrixXd rw2;
    int nmax;
  };
  
  //@{
  //! Grid storage and parameters
  std::vector<RGrid> dens_grid, potl_grid;
  Eigen::VectorXd r_grid;
  double r_grid_del;
  //@}

  //! Cache roots for spherical Bessel functions
  class Roots
  {
  public:

    int l;
    int n;

    Eigen::VectorXd a;

    Roots(int L, int nmax) : l(L), n(nmax) {
      a = sbessjz(l, n);
    }

    ~Roots() {}
  };

  //! Root database isntance
  std::shared_ptr<Roots> p;

  //@{
  //! Density and potential members
  double dens(double r, int n);
  double potl(double r, int n);
  //@}

  //@{
  //! Parameters
  double rmin=0.0, rmax;
  int lmax, nmax, RNUM=2000;
  //@}

  //! Make the density and potential grid
  void make_grid();
  
public:

  //! Constructor
  BiorthBess(double rmax, int lmax, int nmax, int rnum);

  //! Destructor
  virtual ~BiorthBess() {}
	
  //! Potential evaluation
  void get_potl(double r, Eigen::MatrixXd& p);

  //! Density evaluation
  void get_dens(double r, Eigen::MatrixXd& p);

  //! Force evaluation
  void get_dpotl(double r, Eigen::MatrixXd& p);

  //! Potential and density evaluation
  void get_potl_dens(double r, Eigen::MatrixXd& p, Eigen::MatrixXd& d);

  //! Compute the orthogonality of the basis by returning inner
  //! produce matrices
  std::vector<Eigen::MatrixXd> orthoCheck(int knots=40);
};

#endif

