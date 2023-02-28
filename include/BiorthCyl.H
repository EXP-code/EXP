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

//!! BiorthCyl grid class
class BiorthCyl
{

protected:

  YAML::Node conf;

  std::string geometry, forceID;
  
  int mmax, nmax, numr, nfid, mmin, mlim, nmin, nlim, knots, NQDHT;
  double rmin, rmax, scale, acyl;

  bool EVEN_M, verbose;
  
  //@{
  //! Grid parameters
  int cmapR, cmapZ, numx, numy;
  double xmin, xmax, dx, Rtable;
  double ymin, ymax, dy;
  //@}

  //! Storage for basis arrays
  std::vector<std::vector<Eigen::MatrixXd>> dens, pot, rforce, zforce;

  //! Create storage arrays
  void initialize();

  //! Build basis tables
  void create_tables();

  //! Interpolate on grid
  double interp(int m, int n, double R, double z,
		const std::vector<std::vector<Eigen::MatrixXd>>& mat);

  //! Matrix interpolation on grid for coefficient composition
  void interp(double R, double z,
	      const std::vector<std::vector<Eigen::MatrixXd>>& mat,
	      Eigen::MatrixXd& ret);

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

  //! Table scaling factors
  double pfac, ffac, dfac;

public:

  //! Flag for MPI enabled (default: 0=off)
  static int mpi;

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
  double z_to_y(double z);

  //! Map of nondimensional coordinate to vertical coordinate
  double y_to_z(double y);

  //! Jacobian of nondimensional vertical coordinate mapping
  double d_y_to_z(double y);
  
  //! Get potential for dimensionless coord with harmonic order m and radial orer n
  double get_pot(double r, double z, int m, int n)
  { return interp(m, n, r, z, pot) * pfac; }

  //! Get density for dimensionless coord with harmonic order l and radial orer n  
  double get_dens(double r, double z, int m, int n)
  { return interp(m, n, r, z, dens) * dfac; }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  double get_rforce(double r, double z, int m, int n)
  { return interp(m, n, r, z, rforce) * ffac; }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  double get_zforce(double r, double z, int m, int n)
  { return interp(m, n, r, z, zforce) * ffac; }

  //! Get potential for dimensionless coord with harmonic order m and radial orer n
  void get_pot(Eigen::MatrixXd& p, double r, double z)
  { interp(r, z, pot, p); p *= pfac; }

  //! Get density for dimensionless coord with harmonic order l and radial orer n  
  void get_dens(Eigen::MatrixXd& d, double r, double z)
  { interp(r, z, dens, d); d *= dfac; }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  void get_rforce(Eigen::MatrixXd& f, double r, double z)
  { interp(r, z, rforce, f); f *= ffac; }

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  void get_zforce(Eigen::MatrixXd& f, double r, double z)
  { interp(r, z, zforce, f); f *= ffac; }

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

  //@{
  //! Get the current minimum and maximum radii for the expansion
  double getRmin()   { return rmin;   }
  double getRmax()   { return rmax;   }
  double getRtable() { return Rtable; }
  double getRscale() { return scale;  }
  //@}

};

#endif
