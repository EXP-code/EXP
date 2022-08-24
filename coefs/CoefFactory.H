#ifndef _COEF_FACTORY_H
#define _COEF_FACTORY_H

// Needed by member functions for writing parameters and stanzas
#include <highfive/H5File.hpp>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

// The EXP native coefficient classes
#include <CoefsE.H>

/** 
    Abstract class for any time of coefficient database

    Can read and write EXP native and EXP H5 coefficient formats.
    Provides access to coefficients by time. Provides a factory member
    that returns derived class instances.
 */
class Coefs
{

protected:

  //! Blank instance
  Eigen::MatrixXcd mat;

  //! Working array for power
  Eigen::MatrixXd power;

  //! Coefficient type id
  std::string coefType;

  //! Verbose debugging output
  bool verbose;

  //! Read the coefficients
  virtual void readNativeCoefs(const std::string& file) = 0;

  //! Get the YAML config
  virtual std::string getYAML() = 0;

  //! Write parameter attributes
  virtual void WriteH5Params(HighFive::File& file) = 0;

  //! Write coefficient data in H5
  virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count) = 0;

  //! Round time key
  inline double roundTime(double time)
  {
    const double p10 = 1.0e-08;	// 8 signficant figures in coefficient map
    unsigned T = floor(time/p10);
    return p10*T;
  }

public:

  //! Constructor
  Coefs(std::string coefType, bool verbose) :
    coefType(coefType), verbose(verbose) {}

  //! Get coefficient matrix at given time
  virtual Eigen::MatrixXcd& operator()(double time) = 0;

  //! Dump to ascii list
  virtual void dump(int mmin, int mmax, int nmin, int nmax) = 0;

  //! Get list of coefficient times
  virtual std::vector<double> Times() = 0;

  //! Write H5 coefficient file
  virtual void WriteH5Coefs(const std::string& prefix);

  //! Add to an H5 coefficient file
  virtual void ExtendH5Coefs(const std::string& prefix);

  /** Get power for the coefficient DB as a function of harmonic
      index.  Time as rows, harmonics as columns.  This could be
      refined for partial summations. */
  virtual Eigen::MatrixXd& Power() = 0;

  //! The coefficient factory.  Creates a shared pointer to derived
  //! instance of the desired type and returns it.
  static std::shared_ptr<Coefs> factory(const std::string& file);

  //! Compare two collections of stanzas (this is for testing only)
  virtual bool CompareStanzas(std::shared_ptr<Coefs> check) = 0;
};

/** Derived class for spherical coefficients */
class SphCoefs : public Coefs
{
private:
  std::map<double, SphCoefsEPtr> coefs;

protected:

  //! Read the coefficients
  virtual void readNativeCoefs(const std::string& file);

  //! Get the YAML config
  virtual std::string getYAML();

  //! Write parameter attributes
  virtual void WriteH5Params(HighFive::File& file);

  //! Write coefficient data in H5
  virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count);

  /** Get power for the coefficient DB as a function of harmonic
      index.  Time as rows, harmonics as columns.  This could be
      refined for partial summations. */
  virtual Eigen::MatrixXd& Power();


public:

  //! Constructor
  SphCoefs(bool verbose=false) : Coefs("Sphere", verbose) {}

  //! H5 constructor
  SphCoefs(HighFive::File& file, bool verbose=false);

  //! Get coefficient matrix at given time
  virtual Eigen::MatrixXcd& operator()(double time);

  //! Dump to ascii list
  virtual void dump(int lmin, int lmax, int nmin, int nmax);

  //! Get list of coefficient times
  virtual std::vector<double> Times()
  {
    std::vector<double> ret;
    for (auto t : coefs) ret.push_back(t.first);
    return ret;
  }

  //! Compare two collections of stanzas (this is for testing only)
  virtual bool CompareStanzas(std::shared_ptr<Coefs> check);
};

/** Derived class for cylindrical coefficients */
class CylCoefs : public Coefs
{
private:
  std::map<double, CylCoefsEPtr> coefs;
  bool angle;

protected:

  //! Read the coefficients
  virtual void readNativeCoefs(const std::string& file);

  //! Get the YAML config
  virtual std::string getYAML();

  //! Write parameter attributes
  virtual void WriteH5Params(HighFive::File& file);

  //! Write coefficient data in H5
  virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count);

public:

  //! Constructor
  CylCoefs(bool verbose=false) : Coefs("Cylinder", verbose), angle(false) {}

  //! H5 constructor
  CylCoefs(HighFive::File& file, bool verbose=false);

  //! Get coefficient matrix at given time
  virtual Eigen::MatrixXcd& operator()(double time);

  //! Dump to ascii list
  virtual void dump(int mmin, int mmax, int nmin, int nmax);

  //! Get list of coefficient times
  virtual std::vector<double> Times()
  {
    std::vector<double> ret;
    for (auto t : coefs) ret.push_back(t.first);
    return ret;
  }

  /** Get power for the coefficient DB as a function of harmonic
      index.  Time as rows, harmonics as columns.  This could be
      refined for partial summations. */
  virtual Eigen::MatrixXd& Power();

  //! Compare two collections of stanzas (this is for testing only)
  virtual bool CompareStanzas(std::shared_ptr<Coefs> check);
};

using CoefPtr = std::shared_ptr<Coefs>;

#endif
