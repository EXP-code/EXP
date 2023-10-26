#ifndef _BiorthCube_h
#define _BiorthCube_h

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

//!! BiorthCube grid class
class BiorthCube
{

protected:

  YAML::Node conf;

  std::string geometry, forceID;
  
  int nmaxx, nmaxy, nmaxz;
  int nminx, nminy, nminz;
  
  int knots;

public:

  //! Coefficient type
  using coefType = Eigen::Tensor<std::complex<double>, 3>;

  //! Flag for MPI enabled (default: 0=off)
  static int mpi;

  //! Constructor
  BiorthCube(const YAML::Node& conf);

  //! Destructor
  virtual ~BiorthCube() {}

  //! Get potential for dimensionless coord with harmonic order m and radial orer n
  std::complex<double> pot(Eigen::Vector3d x, Eigen::Vector3i n);

  //! Get density for dimensionless coord with harmonic order l and radial orer n  
  std::complex<double> dens(Eigen::Vector3d x, Eigen::Vector3i n);

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  Eigen::Vector3cd force(Eigen::Vector3d x, Eigen::Vector3i n);

  //! Get potential for dimensionless coord with harmonic order m and radial orer n
  std::complex<double> get_pot(const coefType& c, Eigen::Vector3d x);

  //! Get density for dimensionless coord with harmonic order l and radial orer n  
  std::complex<double> get_dens(const coefType& c, Eigen::Vector3d x);

  //! Get radial force for dimensionless coord with harmonic order l and radial orer n
  Eigen::Vector3cd get_force(const coefType& c, Eigen::Vector3d x);

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

  //! For pyEXP
  Eigen::MatrixXcd orthoCheck();

};

#endif
