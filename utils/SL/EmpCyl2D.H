#ifndef _EmpCyl2D_H
#define _EmpCyl2D_H

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <gaussQ.H>

class EmpCyl2D
{
protected:

  double A, rmin, rmax, scale;
  int mmax, nmax, numr, knots;
  std::string model;
  bool cmap, logr;

  std::vector<Eigen::MatrixXd> potl_array, dens_array, dpot_array;
  std::vector<Eigen::MatrixXd> rot_matrix;
  Eigen::VectorXd xgrid;

  double CB_get_potl(int M, int N, double r);
  double CB_get_dens(int M, int N, double r);
  double CB_get_dpot(int M, int N, double r);
  double CB_norm    (int n, int m);
  
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
  
public:

  static std::string cache_name_2d;

  EmpCyl2D(int mmax, int nmax, int knots, int numr,
	   double rmin, double rmax, double A, double scale,
	   bool cmap, bool logr, const std::string type);
  
  void writeBasis(int M, const std::string& file);
  void writeTrans(int M, const std::string& file);
  void orthoCheck(int M, const std::string& file);
};


#endif
