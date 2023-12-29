#ifndef _BiorthCyl_h
#define _BiorthCyl_h

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <Eigen/Eigen>

#include <mpi.h>
#include <localmpi.H>

#include <config_exp.h>		// EXP configuration

#include <massmodel.H>
#include <yaml-cpp/yaml.h>

#if HAVE_LIBCUDA==1
#include <cudaUtil.cuH>
#include <cudaMappingConstants.cuH>
#endif

// For reading and writing cache file
//
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

#if HAVE_LIBCUDA==1
#include <cudaParticle.cuH>
#include <cudaMappingConstants.cuH>
#endif

#include <EmpCyl2d.H>

//!! BiorthCyl grid class
class BiorthCyl
{

protected:

  YAML::Node conf;

  std::string geometry, forceID;
  
  int mmax, nmax, numr, nmaxfid, mmin, mlim, nmin, nlim, knots, NQDHT;
  double rcylmin, rcylmax, scale, acyltbl;

  bool EVEN_M, verbose;
  
  //@{
  //! Grid parameters
  int cmapR, cmapZ, numx, numy;
  double xmin, xmax, dx;
  double ymin, ymax, dy;
  //@}

  //! Storage for basis arrays
  std::vector<std::vector<Eigen::MatrixXd>> dens, pot, rforce, zforce;

  //! The 2d basis instance
  EmpCyl2d emp;

  //! Create storage arrays
  void initialize();

  //! Build basis tables
  void create_tables();

  //! Interpolate on grid
  double interp(int m, int n, double R, double z,
		const std::vector<std::vector<Eigen::MatrixXd>>& mat,
		bool anti_symmetric=false);

  //! Matrix interpolation on grid for coefficient composition
  void interp(double R, double z,
	      const std::vector<std::vector<Eigen::MatrixXd>>& mat,
	      Eigen::MatrixXd& ret, bool anti_symmetric=false);

  //! Default cache file name
  const std::string default_cache = ".biorth_cache";

  //! Active cache file name
  std::string cachename;

  //! Write parameter attributes (needed for derived classes)
  virtual void WriteH5Params(HighFive::File& file);
    
  //! Write the basis arrays
  virtual void WriteH5Arrays(HighFive::Group& arrays);

  //! Read the basis arrays
  virtual void ReadH5Arrays(HighFive::Group& arrays);

  //! Write the HDF5 cache
  virtual void WriteH5Cache();

  //! Read the HDF5 cache
  virtual bool ReadH5Cache();

  //! Flag for MPI enabled (default: 0=off)
  int use_mpi;

public:

  //! Constructor
  BiorthCyl(const YAML::Node& conf);

  //! Destructor
  virtual ~BiorthCyl() {}

  //! Read the cache and report parameters
  static YAML::Node getHeader(const std::string& cachefile);

  //! Map radial coordinate to nondimensional coordinate
  double r_to_xi(double r);

  //! Map nondimensional coordinate to radial coordinate
  double xi_to_r(double x);

  //! Jacobian of nondimensional radial coordinate mapping
  double d_xi_to_r(double x);

  //! Map of vertical coordinate to nondimesional coordinate
  double z_to_yi(double z);

  //! Map of nondimensional coordinate to vertical coordinate
  double yi_to_z(double y);

  //! Jacobian of nondimensional vertical coordinate mapping
  double d_yi_to_z(double y);
  
  //! Get potential for dimensionless coord with harmonic order m and radial orer n
  double get_pot(double r, double z, int m, int n)
  { return interp(m, n, r, z, pot); }

  //! Get density for dimensionless coord with harmonic order l and radial orer n  
  double get_dens(double r, double z, int m, int n)
  { return interp(m, n, r, z, dens); }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  double get_rforce(double r, double z, int m, int n)
  { return interp(m, n, r, z, rforce); }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  double get_zforce(double r, double z, int m, int n)
  { return interp(m, n, r, z, zforce, true); }

  //! Get potential for dimensionless coord with harmonic order m and radial orer n
  void get_pot(Eigen::MatrixXd& p, double r, double z)
  { interp(r, z, pot, p); }

  //! Get density for dimensionless coord with harmonic order l and radial orer n  
  void get_dens(Eigen::MatrixXd& d, double r, double z)
  { interp(r, z, dens, d); }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  void get_rforce(Eigen::MatrixXd& f, double r, double z)
  { interp(r, z, rforce, f); }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  void get_zforce(Eigen::MatrixXd& f, double r, double z)
  { interp(r, z, zforce, f, true); }

  //! Read and print the cache and return the header parameters as a
  //! map/dictionary
  static std::map<std::string, std::string>
  cacheInfo(const std::string& cachefile, bool verbose=true);

#if HAVE_LIBCUDA==1
  void initialize_cuda(std::vector<cudaArray_t>& cuArray,
		       thrust::host_vector<cudaTextureObject_t>& tex);

  virtual cudaMappingConstants getCudaMappingConstants()
  {
    cudaMappingConstants ret;

    ret.rscale = scale;
    ret.hscale = scale;
    ret.xmin   = xmin;
    ret.xmax   = xmax;
    ret.ymin   = ymin;
    ret.ymax   = ymax;
    ret.numr   = numr;
    ret.numx   = numx;
    ret.numy   = numy;
    ret.dxi    = dx;
    ret.dyi    = dy;
    ret.cmapR  = cmapR;
    ret.cmapZ  = cmapZ;

    return ret;
  }
#endif

  //! Evaluate all orders in matrices; for n-body
  void get_pot(Eigen::MatrixXd& Vc, Eigen::MatrixXd& Vs, double r, double z);

  //! Get the table bounds
  double getRtable() { return rcylmax*scale; }

  //! For pyEXP
  std::vector<Eigen::MatrixXd> orthoCheck();

};

#endif
