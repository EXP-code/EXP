#ifndef _COEFS_H
#define _COEFS_H

#include <memory>
#include <Eigen/Eigen>

#include <yaml-cpp/yaml.h>	// YAML support

#include "coef.H"

struct CylCoefs
{
  std::shared_ptr<char[]> buf;
  double time;
  int nmax, mmax;
  std::vector< std::vector<double> > cos_c, sin_c;

  bool read(std::istream& in, bool verbose=false);
};

typedef std::shared_ptr<CylCoefs> CylCoefsPtr;

struct SphCoefs
{
  std::shared_ptr<char[]> buf;
  SphCoefHeader header;

  typedef
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
  EigenColMajor;

  EigenColMajor coefs;

  bool read(std::istream& in, bool exp_type);
};

typedef std::shared_ptr<SphCoefs> SphCoefsPtr;

#endif
