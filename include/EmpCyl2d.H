#ifndef _EmpCyl2d_H
#define _EmpCyl2d_H

#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <tuple>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>

#include <yaml-cpp/yaml.h>	// Parameters parsing
#include "config_exp.h"		// For config macros

#include "gaussQ.H"
#include "EXPmath.H"
#include "localmpi.H"

// For reading and writing cache file
//
#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>

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

    std::string id;
    virtual void parse(const YAML::Node&) = 0;

  public:

    virtual double pot (double r) = 0;
    virtual double dpot(double r) = 0;
    virtual double dens(double r) = 0;

    virtual std::string ID() { return id; }
  };

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


  class ExponCyl : public ModelCyl
  {
    
  private:
    
    double acyl, sigma0;

    // Assign values from YAML
    //
    void parse(const YAML::Node& conf)
    {
      try {
	if (conf["acyl"]) 
	  acyl = conf["acyl"].as<double>();
	else
	  acyl = 1.0;
      }
      catch (YAML::Exception & error) {
	if (myid==0)
	  std::cout << "Error parsing parameters in EmpCyl2d::ExponCyl: "
		    << error.what() << std::endl
		    << std::string(60, '-') << std::endl
		    << "Config node"        << std::endl
		    << std::string(60, '-') << std::endl
		    << conf                 << std::endl
		    << std::string(60, '-') << std::endl;
	throw std::runtime_error("EmpCyl2d::ExponCyl: error parsing YAML config");
      }
    }
    
  public:

    ExponCyl(const YAML::Node& par)
    {
      parse(par);
      sigma0 = 0.5/(M_PI*acyl*acyl);
      id = "expon";
    }
    
    double pot(double r) {
      double y = 0.5 * r / acyl;
      return -M_PI*sigma0*r *
	(EXPmath::cyl_bessel_i(0, y)*EXPmath::cyl_bessel_k(1, y) -
	 EXPmath::cyl_bessel_i(1, y)*EXPmath::cyl_bessel_k(0, y));
    }

    double dpot(double r) {
      double y = 0.5 * r / acyl;
      return 2.0*M_PI*sigma0*y*
	(EXPmath::cyl_bessel_i(0, y)*EXPmath::cyl_bessel_k(0, y) -
	 EXPmath::cyl_bessel_i(1, y)*EXPmath::cyl_bessel_k(1, y));
    }
    
    double dens(double r) {
      return sigma0*exp(-r/acyl);
    }
    
  };

  class KuzminCyl : public ModelCyl
  {
  private:

    double acyl;

    // Assign values from YAML
    //
    void parse(const YAML::Node& conf)
    {
      try {
	if (conf["acyl"]) 
	  acyl = conf["acyl"].as<double>();
	else
	  acyl = 1.0;
      }
      catch (YAML::Exception & error) {
	if (myid==0)
	  std::cout << "Error parsing parameters in EmpCyl2d::KuzminCyl: "
		    << error.what() << std::endl
		    << std::string(60, '-') << std::endl
		    << "Config node"        << std::endl
		    << std::string(60, '-') << std::endl
		    << conf                 << std::endl
		    << std::string(60, '-') << std::endl;
	throw std::runtime_error("EmpCyl2d::KuzminCyl: error parsing YAML config");
      }
    }
    
  public:
  
    KuzminCyl(const YAML::Node& par)
    {
      parse(par);
      id = "kuzmin";
    }

    double pot(double R) {
      double a2 = acyl * acyl;
      return -1.0/sqrt(R*R + a2);
    }
    
    double dpot(double R) {
      double a2 = acyl * acyl;
      return R/pow(R*R + a2, 1.5);
    }
    
    double dens(double R) {
      double a2 = acyl * acyl;
      return 4.0*M_PI*acyl/pow(R*R + a2, 1.5)/(2.0*M_PI);
      //     ^
      //     |
      // This 4pi from Poisson's eqn
    }
    
  };


  class MestelCyl : public ModelCyl
  {
  protected:

    double vrot, rot;

    // Assign values from YAML
    //
    void parse(const YAML::Node& conf)
    {
      try {
	if (conf["vrot"]) 
	  vrot = conf["vrot"].as<double>();
	else
	  vrot = 1.0;
      }
      catch (YAML::Exception & error) {
	if (myid==0)
	  std::cout << "Error parsing parameters in EmpCyl2d::MestelCyl: "
		    << error.what() << std::endl
		    << std::string(60, '-') << std::endl
		    << "Config node"        << std::endl
		    << std::string(60, '-') << std::endl
		    << conf                 << std::endl
		    << std::string(60, '-') << std::endl;
	throw std::runtime_error("EmpCyl2d::MestelCyl: error parsing YAML config");
      }
    }
    
  public:
  
    MestelCyl(const YAML::Node& par)
    {
      parse(par);
      rot = vrot*vrot;
      id = "mestel";
    }

    virtual double pot(double R) {
      if (R>0.0) return rot*log(R);
      else throw std::runtime_error("MestelCyl::pot: R<=0");
    }
    
    virtual double dpot(double R) {
      if (R>0.0) return rot/R;
      else throw std::runtime_error("MestelCyl::dpot: R<=0");
    }

    virtual double dens(double R) {
      if (R>0.0) return rot/(2.0*M_PI*R);
      else throw std::runtime_error("MestelCyl::dens: R<=0");
    }
  };
  

  class ZangCyl : public MestelCyl
  {
  
  private:
    //! Parameters
    double mu, nu, ri, ro;

    //! Softening factor (not currently used)
    double asoft = 1.0e-8;

    //! Ignore inner cut-off for N<0.05
    bool Inner = true;

    //! Taper factors
    double Tifac, Tofac;

    //! Inner taper function
    double Tinner(double Jp)
    {
      double fac = pow(Jp, nu);
      return fac/(Tifac + fac);
    }

    //! Outer taper function
    double Touter(double Jp)
    {
      return 1.0/(1.0 + pow(Jp/Tofac, mu));
    }
    
    //! Deriv of inner taper function
    double dTinner(double Jp)
    {
      double fac  = pow(Jp, nu);
      double fac2 = Tifac + fac;
      return nu*fac/Jp/(fac2*fac2);
    }
    
    //! Deriv of outer taper function
    double dTouter(double Jp)
    {
      double fac = pow(Jp/Tofac, mu);
      double fac2 = 1.0 + fac;
      return -mu*fac/Jp/(fac2*fac2);
    }
    
  protected:
    
    //! Assign values from YAML
    void parse(const YAML::Node& conf)
    {
      try {
	if (conf["Ninner"]) 
	  nu = conf["Ninner"].as<double>();
	else
	  nu = 2.0;

	if (conf["Mouter"]) 
	  mu = conf["Mouter"].as<double>();
	else
	  mu = 2.0;
	
	if (conf["Ri"]) 
	  ri = conf["Ri"].as<double>();
	else
	  ri = 1.0;
	
	if (conf["Ro"]) 
	  ro = conf["Ro"].as<double>();
	else
	  ro = 10.0;
      }
      catch (YAML::Exception & error) {
	if (myid==0)
	  std::cout << "Error parsing parameters in EmpCyl2d::ZangCyl: "
		    << error.what() << std::endl
		    << std::string(60, '-') << std::endl
		    << "Config node"        << std::endl
		    << std::string(60, '-') << std::endl
		    << conf                 << std::endl
		    << std::string(60, '-') << std::endl;
	throw std::runtime_error("EmpCyl2d::ZangCyl: error parsing YAML config");
      }
    }

  public:
  
    //! Constructor
    ZangCyl(const YAML::Node& par) : MestelCyl(par)
    {
      // Parse the YAML
      parse(par);
      
      // Assign the id
      id = "zang";
      
      // Cache taper factors
      Tifac = pow(ri*vrot, nu);
      Tofac = ro*vrot;
      
      if (nu<0.05) {
	// Exponent is now for mapping only
	Inner = false;
      }
    }
    
    //! Surface density
    double dens(double R)
    {
      double ret = MestelCyl::dens(R) * Touter(R);
      if (Inner) ret *= Tinner(R);
      return ret;
    }
    
  };


protected:

  //! Contains parameter database
  YAML::Node Params;

  double rmin, rmax, scale;
  double xmin, xmax, dxi;
  int mmax, nmaxfid, numr, knots, nmax;
  std::string model, biorth;
  bool cmap, logr, basis_test, configured;

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
    createBasis(int mmax, int nmaxfid, double rmax, const std::string& type);

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

  std::shared_ptr<ModelCyl> createModel();

  //! Mapping instance
  Mapping map;

  void create_tables();
  bool ReadH5Cache();
  void WriteH5Cache();

  //! Cache versioning
  inline static const std::string Version = "1.0";

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
  EmpCyl2d() : configured(false) {}

  //! Constructor
  EmpCyl2d(int mmax, int nmaxfid, int nmax, int knots, int numr,
	   double rmin, double rmax, double scale, bool cmap, bool logr,
	   const YAML::Node& P,
	   const std::string biorth, const std::string cache="");
  
  //! Constructor with user-supplied target model
  EmpCyl2d(int mmax, int nmaxfid, int nmax, int knots, int numr,
	   double rmin, double rmax, double scale, bool cmap, bool logr,
	   std::shared_ptr<ModelCyl> disk,
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

  //! Check for configuration
  bool operator()() { return configured; }

  //! Get coordinate mapping object
  Mapping getMapping() { return map; }

  //! Get model name
  const std::string getModelName() { return model; }

  //! Get biorthogonal function base name
  const std::string getBiorthName() { return biorth; }

};


#endif
