#ifndef _COEFSE_H
#define _COEFSE_H

#include <memory>
#include <Eigen/Eigen>

//! This is a general struct for any dimensionality
struct CoefsE
{
  std::shared_ptr<char[]> buf;
  double time;

  typedef
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
  EigenColMajor;

  EigenColMajor coefs;

  virtual bool read(std::istream& in, bool exp_type) = 0;
};

//! Specialization of CoefsE for spheres
class SphCoefsE : public CoefsE
{
public:
  //! Angular and radial dimension
  int lmax, nmax;

  //! Radial scale factor
  double scale;

  //! Is basis unit normed?
  bool normed;

  //! Read native EXP coefficients
  bool read(std::istream& in, bool exp_type);
};

//! Specialization of CoefsE for cylinders
class CylCoefsE : public CoefsE
{
public:
  //! Angular and radial dimension
  int mmax, nmax;

  //! Read native EXP coefficients
  bool read(std::istream& in, bool exp_type);
};

typedef std::shared_ptr<CoefsE   > CoefsEPtr;
typedef std::shared_ptr<SphCoefsE> SphCoefsEPtr;
typedef std::shared_ptr<CylCoefsE> CylCoefsEPtr;

#endif
