#ifndef _BiorthCube_h
#define _BiorthCube_h

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <Eigen/Eigen>

#include <mpi.h>
#include "localmpi.H"

#include "config_exp.h"		// EXP configuration

#include "massmodel.H"
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

#include "EmpCyl2d.H"

//!! BiorthCube grid class
class BiorthCube
{

protected:

  YAML::Node conf;

  std::string geometry, forceID;
  
  Eigen::Vector3i nmin, nmax;
  
  int knots;

  bool verbose;

  static constexpr double nfac = 16.0*M_PI*M_PI*M_PI;
  static constexpr double dfac = 2.0*M_PI;
  static constexpr std::complex<double> kfac{0.0, dfac};

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

  //! For pyEXP
  Eigen::MatrixXcd orthoCheck();

};

#endif
