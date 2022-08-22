#ifndef _H5COEFS_H
#define _H5COEFS5_H

#include <highfive/H5File.hpp>

#include <Coefs.H>

class CoefClient;

/** 
    Abstract class for any time of coefficient database

    Can read and write EXP native and EXP H5 coefficient formats
 */
class Coefs
{
  friend class CoefClient;

protected:

  //! Blank instance
  Eigen::MatrixXcd blank;

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
    unsigned T = floor(time*1.0e8);
    return 1.0e-08*T;
  }

public:

  //! Constructor
  Coefs(std::string coefType, bool verbose) : coefType(coefType), verbose(verbose) {}

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

};

/** Derived class for spherical coefficients */
class SphH5 : public Coefs
{
private:
  std::map<double, SphCoefsPtr> coefs;

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
  SphH5(bool verbose=false) : Coefs("Sphere", verbose) {}

  //! H5 constructor
  SphH5(HighFive::File& file, bool verbose=false);

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
class CylH5 : public Coefs
{
private:
  std::map<double, CylCoefsPtr> coefs;
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
  CylH5(bool verbose=false) : Coefs("Cylinder", verbose), angle(false) {}

  //! H5 constructor
  CylH5(HighFive::File& file, bool verbose=false);

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

};

/** Factory, wrapper class for Coef */
class CoefClient
{
private:
  
  //! Factory instance
  std::shared_ptr<Coefs> coefs;

public:

  //! Constructor
  CoefClient(const std::string& file);

  //! Write to H5
  void WriteH5Coefs(const std::string& file)  { coefs->WriteH5Coefs(file);  }

  //! Extend H5file
  void ExtendH5Coefs(const std::string& file) { coefs->ExtendH5Coefs(file); }

  //! Instance
  Coefs* get() { return coefs.get(); }

  //! Get list of times
  std::vector<double> Times() { return coefs->Times(); }

  //! Compare two coefficient collections (this is for testing only)
  bool CompareStanzas(CoefClient& check)
  { return coefs->CompareStanzas(check.get()); }

  //! Get coefficient matrix at given time
  Eigen::MatrixXcd& operator()(double time) { return (*coefs)(time); }


};

#endif
