#ifndef _COEF_FACTORY_H
#define _COEF_FACTORY_H

// Needed by member functions for writing parameters and stanzas
#include <highfive/H5File.hpp>

#include <Eigen/Eigen>	      // Store coefficient matrices and 2d grids
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids?

// The EXP native coefficient classes
#include <CoefStruct.H>

namespace Coefs
{ 
  /** 
      Abstract class for any type of coefficient database
      
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
    
    //! Get the YAML config for the basis (to be implemented by EXP)
    virtual std::string getYAML() = 0;
    
    //! Write parameter attributes (needed for derived classes)
    virtual void WriteH5Params(HighFive::File& file) = 0;
    
    //! Write coefficient data in H5
    virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count) = 0;
    
    //! Round time key to emulated fixed-point arithmetic
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
    
    //! Get coefficient structure at a given time
    virtual std::shared_ptr<CoefStruct> getCoefStruct(double time) = 0;

    //! Dump to ascii list (for testing only)
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
    
    /** The coefficient factory members

	Creates a shared pointer to derived instance of the desired
	type and returns it.  This uses a coefficient file to
	construct the coefficient database
    */
    static std::shared_ptr<Coefs> factory(const std::string& file);
    
    /** Joint factory and update member

	If Coefs pointer is null:
	* Creates a coefficient container based on the initial set of
	  coefficients and adds them to the internal container

	If Coefs pointer is not null:
	* Use add() member to add new coefficients to the exising
	  container

	This factory would be used with
	Basis::Basis to derive coefficients from snapshots.
    */
    static std::shared_ptr<Coefs>
    addcoef(std::shared_ptr<Coefs> coefs, CoefStrPtr coef);
    
    //! Add a coefficient structure to the container
    virtual void add(CoefStrPtr coef) = 0;
    
    //! Compare two collections of stanzas (this is for testing only)
    virtual bool CompareStanzas(std::shared_ptr<Coefs> check) = 0;
  };
  
  /** Derived class for spherical coefficients */
  class SphCoefs : public Coefs
  {
  private:
    std::map<double, SphStrPtr> coefs;
    
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
    
    //! Add a coefficient structure to the container
    virtual void add(CoefStrPtr coef);

    //! Get coefficient matrix at given time
    virtual Eigen::MatrixXcd& operator()(double time);
    
    //! Get coefficient structure at a given time
    virtual std::shared_ptr<CoefStruct> getCoefStruct(double time)
    { return coefs[time]; }

    //! Dump to ascii list for testing
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
    std::map<double, CylStrPtr> coefs;
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
    
    //! Add a coefficient structure to the container
    virtual void add(CoefStrPtr coef);

    //! Get coefficient matrix at given time
    virtual Eigen::MatrixXcd& operator()(double time);
    
    //! Get coefficient structure at a given time
    virtual std::shared_ptr<CoefStruct> getCoefStruct(double time)
    { return coefs[time]; }


    //! Dump to ascii list for testing
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
  
  using CoefsPtr = std::shared_ptr<Coefs>;
}
// END namespace Coefs
  
#endif
  
