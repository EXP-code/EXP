#ifndef _COEFSTOH5_H
#define _COEFSTOH5_H

#include <highfive/H5File.hpp>

#include <Coefs.H>

class CoefClient;

class Coefs
{
  friend class CoefClient;

protected:

  //! Coefficient type id
  std::string coefType;

  //! Verbose debugging output
  bool verbose;

  //! Read the coefficients
  virtual void readNativeCoefs(const std::string& file) = 0;

  //! Get the YAML config
  virtual std::string getYAML() = 0;

  //! Write coefficient data in H5
  virtual void WriteH5Times(HighFive::Group& group) = 0;

public:

  //! Constructor
  Coefs(std::string coefType, bool verbose) : coefType(coefType), verbose(verbose) {}

  //! Dump to ascii list
  virtual void dump(int mmin, int mmax, int nmin, int nmax) = 0;

  //! Get list of coefficient times
  virtual std::vector<double> Times() = 0;

  //! Write H5 coefficient file
  virtual void WriteH5Coefs(const std::string& prefix);

};

class SphH5 : public Coefs
{
private:
  std::map<double, SphCoefsPtr> coefs;

protected:

  //! Read the coefficients
  virtual void readNativeCoefs(const std::string& file);

  //! Get the YAML config
  virtual std::string getYAML();

  //! Write coefficient data in H5
  virtual void WriteH5Times(HighFive::Group& group);

public:

  //! Constructor
  SphH5(bool verbose=true) : Coefs("Sphere", verbose) {}

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

  //! Write coefficient data in H5
  virtual void WriteH5Times(HighFive::Group& group);

public:

  //! Constructor
  CylH5(bool verbose=true) : Coefs("Cylinder", verbose), angle(false) {}

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

class CoefClient
{
private:
  
  //! Factory instance
  std::shared_ptr<Coefs> coefs;

public:

  //! Constructor
  CoefClient(const std::string& file);

  //! Write to H5
  void WriteH5Coefs(const std::string& file) { coefs->WriteH5Coefs(file); }

  //! Instance
  Coefs* get() { return coefs.get(); }

  //! Get list of times
  std::vector<double> Times() { return coefs->Times(); }

};

#endif
