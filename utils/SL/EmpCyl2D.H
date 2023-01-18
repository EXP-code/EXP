#ifndef _EmpCyl2D_H
#define _EmpCyl2D_H

#include <algorithm>
#include <string>
#include <vector>
#include <tuple>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <gaussQ.H>

/**
   A class that implements most of the members for an Exp force routine

   With a bit of work, this could be wrapped into a Basis class.

   Things necessary to do this:

   1. Add MPI support.  This could be fairly minimal since the
      overhead in this class if very low.  At the very least, only the
      root should be allowed to write a cache file.

   2. Add Cuda support.  This guts could be straightforward copied
      from the 1-d SLGridSph array and interpolation support.
 */
class EmpCyl2D
{
protected:

  double A, rmin, rmax, scale;
  int mmax, nmax, numr, knots;
  std::string model, biorth;
  bool cmap, logr;

  std::vector<Eigen::MatrixXd> potl_array, dens_array, dpot_array;
  std::vector<Eigen::MatrixXd> rot_matrix;
  std::vector<double> xgrid;

  //! Two-d basis models
  class Basis2d
  {
  public:

    static std::shared_ptr<Basis2d>
    createBasis(int mmax, int nmax, const std::string& type);

    virtual double potl(int M, int N, double r) = 0;
    virtual double dens(int M, int N, double r) = 0;
    virtual double dpot(int M, int N, double r) = 0;
    virtual double norm(int N, int M) = 0;
  };

  class CluttonBrock;
  class Bessel;

  std::shared_ptr<Basis2d> basis;

  //! A two-dimensional disk model for computing the EOF
  class ModelCyl
  {
  protected:

    double A;
    std::string id;

  public:

    static std::shared_ptr<ModelCyl>
    createModel(const std::string type, double A);

    virtual double pot (double r) = 0;
    virtual double dpot(double r) = 0;
    virtual double dens(double r) = 0;
  };
  
  class ExponCyl;
  class MestelCyl;
  class KuzminCyl;

  //! Map the radius
  class Mapping
  {
  protected:
    bool cmap;
    double scale;
    
  public:
    
    Mapping(double scale, bool cmap) : scale(scale), cmap(cmap) {}

    double r_to_xi(double r);
    double xi_to_r(double xi);
    double d_xi_to_r(double xi);
  };

  void create_tables();
  bool read_cached_tables();
  void write_cached_tables();

  //! Basis magic number
  inline static const unsigned int hmagic = 0xc0a57a1;
  
  //@{
  //! Interpolation helpers
  std::tuple<int, int, double, double> linear_interp(double r);
  void checkMN(int& M, int& N, const std::string& member);
  //@}
  
public:

  //! The cache file name
  static std::string cache_name_2d;

  //! Constructor
  EmpCyl2D(int mmax, int nmax, int knots, int numr,
	   double rmin, double rmax, double A, double scale,
	   bool cmap, bool logr, const std::string type,
	   const std::string biorth);
  
  //! Linear interpolate on grid to get potential value
  double get_potl(double r, int M, int N);

  //! Linear interpolate on grid to get density value
  double get_dens(double r, int M, int N);

  //! Linear interpolate on grid to get potential gradient (-force) value
  double get_dpot(double r, int M, int N);

  //! Write the basis for a partcular azimuthal M value to a file
  void writeBasis(int M, const std::string& file);

  //! Write the eigenvector rotation matrix for a partcular azimuthal
  //! M value to a file
  void writeTrans(int M, const std::string& file);

  //! Write the inner-product matrix for the EOF basis for a
  //! particular azimuthal M value to a file
  void orthoCheck(int M, const std::string& file);
};


#endif
