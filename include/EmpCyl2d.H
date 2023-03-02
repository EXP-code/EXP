#ifndef _EmpCyl2d_H
#define _EmpCyl2d_H

#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <tuple>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <gaussQ.H>

#include <config_exp.h>		// For config macros

// For reading and writing cache file
//
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

/**
   A class that implements most of the members for an Exp force routine
*/
class EmpCyl2d
{
public:
  //! A two-dimensional disk model for computing the EOF
  class ModelCyl
  {
  protected:

    double acyl;
    std::string id;

  public:

    virtual double pot (double r) = 0;
    virtual double dpot(double r) = 0;
    virtual double dens(double r) = 0;

    virtual std::string ID() { return id; }
  };

protected:

  double acyl, rmin, rmax, scale;
  double xmin, xmax, dxi;
  int mmax, nmax, numr, knots;
  std::string model, biorth;
  bool cmap, logr, basis_test;

  std::vector<Eigen::MatrixXd> potl_array, dens_array, dpot_array;
  std::vector<Eigen::MatrixXd> rot_matrix;
  std::vector<double> xgrid;

  //! Two-d basis models
  class Basis2d
  {
  public:

    //! Enable orthgonality check
    static bool debug;

    //! Basis creation factory
    static std::shared_ptr<Basis2d>
    createBasis(int mmax, int nmax, double rmax, const std::string& type);

    //@{
    //! Required basis-function members
    virtual double potl(int M, int N, double r) = 0;
    virtual double dens(int M, int N, double r) = 0;
    virtual double dpot(int M, int N, double r) = 0;
    virtual double norm(int N, int M) = 0;
    //@}
  };

  class CluttonBrock;
  class Bessel;

  std::shared_ptr<ModelCyl> disk;
  std::shared_ptr<Basis2d>  basis;

  static std::shared_ptr<ModelCyl>
  createModel(const std::string type, double A);

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
    
    //! Null constructor (for copy construct)
    Mapping() {}

    //! Main constructor
    Mapping(double scale, bool cmap) : scale(scale), cmap(cmap) {}

    //! From semi-infinite radius to [-1, 1]
    double r_to_xi(double r);

    //! From [-1, 1] to semi-infinite radius
    double xi_to_r(double xi);

    //! Jacobian
    double d_xi_to_r(double xi);
  };

  //! Mapping instance
  Mapping map;

  void create_tables();
  bool ReadH5Cache();
  void WriteH5Cache();

  //! Basis magic number
  inline static const unsigned int hmagic = 0xc0a57a1;
  
  //@{
  //! Interpolation helpers
  std::tuple<int, int, double, double> linear_interp(double r);
  void checkMN(int& M, int& N, const std::string& member);
  //@}
  
  //! Default cache name
  static const std::string default_cache_name;

  //! The cache file name
  std::string cache_name_2d;

public:

  //! Null constructor (for copy construct)
  EmpCyl2d() {}

  //! Constructor
  EmpCyl2d(int mmax, int nmax, int knots, int numr,
	   double rmin, double rmax, double A, double scale,
	   bool cmap, bool logr, const std::string model,
	   const std::string biorth, const std::string cache="");
  
  //! Constructor with user-supplied target model
  EmpCyl2d(int mmax, int nmax, int knots, int numr,
	   double rmin, double rmax, double A, double scale,
	   bool cmap, bool logr, std::shared_ptr<ModelCyl> disk,
	   const std::string biorth, const std::string cache="");
  
  //! Use underlying basis for testing?  Default: false.
  void basisTest(bool which) { basis_test = which; }

  //! Linear interpolate on grid to get potential value
  double get_potl(double r, int M, int N);

  //! Linear interpolate on grid to get density value
  double get_dens(double r, int M, int N);

  //! Linear interpolate on grid to get potential gradient (-force) value
  double get_dpot(double r, int M, int N);

  //@{
  //! Matrix evaluation members
  void get_pot  (Eigen::MatrixXd& mat, double r);
  void get_dens (Eigen::MatrixXd& mat, double r);
  void get_force(Eigen::MatrixXd& mat, double r);
  //@}

  //! Write the basis for a partcular azimuthal M value to a file
  void writeBasis(int M, const std::string& file);

  //! Write the eigenvector rotation matrix for a partcular azimuthal
  //! M value to a file
  void writeTrans(int M, const std::string& file);

  //! Write the inner-product matrix for the EOF basis for a
  //! particular azimuthal M value to a file
  void orthoCheck(int M, const std::string& file);

  //! Return an array of inner-product matrices for the EOF basis
  std::vector<Eigen::MatrixXd> orthoCheck();

  //! Check the coefficients of the EOF basis against the zeroth-order
  //! model
  void checkCoefs();

  //! Access the minimum radius
  double getRmin() { return rmin; }

  //! Access the maximum radius
  double getRmax() { return rmax; }
};


#endif
