#ifndef _COEF_FACTORY_H
#define _COEF_FACTORY_H

// Needed by member functions for writing parameters and stanzas
#include <highfive/H5File.hpp>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

// The EXP native coefficient classes
#include <CoefsE.H>

// The wrapper class
class CoefFactory;

/** 
    Abstract class for any time of coefficient database

    Can read and write EXP native and EXP H5 coefficient formats.  Provides
    access to coefficients by time.  
 */
class Coefs
{
  friend class CoefFactory;

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

  //! Compare two collections of stanzas (this is for testing only)
  virtual bool CompareStanzas(Coefs* check) = 0;

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
  virtual Eigen::MatrixXd& getPower() = 0;

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

  //! Compare two collections of stanzas (this is for testing only)
  virtual bool CompareStanzas(Coefs* check);

  /** Get power for the coefficient DB as a function of harmonic
      index.  Time as rows, harmonics as columns.  This could be
      refined for partial summations. */
  virtual Eigen::MatrixXd& getPower();


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

  //! Compare two collections of stanzas (this is for testing only)
  virtual bool CompareStanzas(Coefs* check);
  
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
  virtual Eigen::MatrixXd& getPower();

};

/** 
    Factory and wrapper class for Coef

    Provides the API of Coef while hiding the details of the instance.
    The type of coefficient is deduced from the coefficient file.
 */
class CoefFactory
{
private:
  
  //! Factory instance
  std::shared_ptr<Coefs> coefs;

public:

  //! Constructor
  //! @param file can be either a native or HDF5 coefficient file
  CoefFactory(const std::string& file);

  //! Write a new HDF5 coefficient file
  void WriteH5Coefs(const std::string& file)  { coefs->WriteH5Coefs(file);  }

  //! Extend a current HDF5 coefficient file with additional stanzas
  //! from the native file This is currently only a test for the
  //! 'extend' concept and probably needs modification for production.
  void ExtendH5Coefs(const std::string& file) { coefs->ExtendH5Coefs(file); }

  //! Get the underlying Coefs instance
  Coefs* get() { return coefs.get(); }

  //! Get list of times
  std::vector<double> Times() { return coefs->Times(); }

  //! Get total power as a function of harmonic index (type specific)
  Eigen::MatrixXd& Power() { return coefs->getPower(); }

  //! Compare two coefficient collections (this is for testing only)
  bool CompareStanzas(CoefFactory& check)
  { return coefs->CompareStanzas(check.get()); }

  //! Get coefficient matrix at given time
  Eigen::MatrixXcd& operator()(double time) { return (*coefs)(time); }

};

#endif
