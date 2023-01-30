#ifndef _COEFFICIENTS_H
#define _COEFFICIENTS_H

#include <tuple>

// Needed by member functions for writing parameters and stanzas
#include <highfive/H5File.hpp>

#include <Eigen/Eigen>	      // Store coefficient matrices and 2d grids
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids?

// The EXP native coefficient classes
#include <CoefStruct.H>

namespace CoefClasses
{ 
  //! An index key
  using Key = std::vector<unsigned>;
  
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
    
    //! Coefficient geometry
    std::string geometry;
    
    //! Coefficient name (mnemonic)
    std::string name;

    //! Verbose debugging output
    bool verbose;
    
    //! Time vector
    std::vector<double> times;

    //! Read the coefficients
    virtual void readNativeCoefs(const std::string& file,
				 int stride, double tmin, double tmax) = 0;
    
    //! Get the YAML config for the basis (to be implemented by EXP)
    virtual std::string getYAML() = 0;
    
    //! Write parameter attributes (needed for derived classes)
    virtual void WriteH5Params(HighFive::File& file) = 0;
    
    //! Write coefficient data in H5
    virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count) = 0;
    
    //! Round time key to emulated fixed-point arithmetic
    inline double roundTime(double time)
    {
      // Eight decimal places should be enough here...
      const double multiplier = 1.0e+08;  // std::pow(10.0, 8);
      return std::floor(time * multiplier + 0.5) / multiplier;
    }
    
    //! Copy all of the base-class data
    void copyfields(std::shared_ptr<Coefs> p);

    //! Time offset for interpolation
    double deltaT;

  public:
    
    //! Constructor
    Coefs(std::string geometry, bool verbose) :
      geometry(geometry), verbose(verbose), deltaT(0.01) {}

    //! Copy constructor
    Coefs(Coefs& p)
    {
      mat      = p.mat;
      power    = p.power;
      geometry = p.geometry;
      name     = p.name;
      times    = p.times;
      deltaT   = p.deltaT;
    }
    
    //! Get coefficient matrix at given time
    virtual Eigen::MatrixXcd& operator()(double time) = 0;
    
    //! Set coefficient matrix at given time to the provided matrix
    virtual void setMatrix(double time, const Eigen::MatrixXcd& mat) = 0;

    //! Interpolate coefficient matrix at given time
    virtual std::tuple<Eigen::MatrixXcd&, bool> interpolate(double time);
    
    //! Get coefficient structure at a given time
    virtual std::shared_ptr<CoefStruct> getCoefStruct(double time) = 0;

    //! Dump to ascii list (for testing only)
    virtual void dump(int mmin, int mmax, int nmin, int nmax) = 0;
    
    //! Get list of coefficient times
    virtual std::vector<double> Times() { return times; }
    
    //! Write H5 coefficient file
    virtual void WriteH5Coefs(const std::string& prefix);
    
    //! Add to an H5 coefficient file
    virtual void ExtendH5Coefs(const std::string& prefix);
    
    /** Get power for the coefficient DB as a function of harmonic
	index.  Time as rows, harmonics as columns.

	@param min is the minimum radial order
	@param max is the maximum radial order
    */
    virtual Eigen::MatrixXd& Power
    (int min=0, int max=std::numeric_limits<int>::max()) = 0;
    
    //! Make keys for the remaining indices in a subspace
    virtual std::vector<Key> makeKeys(Key k) = 0;

    /** The coefficient factory members

	Creates a shared pointer to derived instance of the desired
	type and returns it.  This uses a coefficient file to
	construct the coefficient database
    */
    static std::shared_ptr<Coefs> factory
    (const std::string& file, int stride=1,
     double tmin=-std::numeric_limits<double>::max(),
     double tmax= std::numeric_limits<double>::max() );
    
    //! Make Coefs instance if it doesn't yet exist
    static std::shared_ptr<Coefs> makecoefs(CoefStrPtr coef, std::string name="");

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
    
    //! Clear coefficient container
    virtual void clear() = 0;

    //! Add a coefficient structure to the container
    virtual void add(CoefStrPtr coef) = 0;
    
    //! Compare two collections of stanzas (this is for testing only)
    virtual bool CompareStanzas(std::shared_ptr<Coefs> check) = 0;

    //! Get coefficient type
    std::string getGeometry() { return geometry; }

    //! Get mnemonic name
    std::string getName() { return name; }

    //! Set mnemonic name
    void setName(const std::string& newname) { name = newname; }

    //! Copy all of the data; make new instances of shared pointer
    //! objects
    virtual std::shared_ptr<Coefs> deepcopy() = 0;

    //! Zero the existing data
    virtual void zerodata() = 0;

    //! Set maximum grid interpolation offset
    void setDeltaT(double dT) { deltaT = dT; }

    class CoefsError : public std::runtime_error
    {
    public:
      CoefsError(const std::string& message) : std::runtime_error(message) {}
    };

  };
  
  /** Derived class for spherical coefficients */
  class SphCoefs : public Coefs
  {
  protected:
    
    //! Parameters
    int Lmax, Nmax;

    //! The coefficient DB
    std::map<double, SphStrPtr> coefs;

    //! Read the coefficients
    virtual void readNativeCoefs(const std::string& file,
				 int stride, double tmin, double tmax);
    
    //! Get the YAML config
    virtual std::string getYAML();
    
    //! Write parameter attributes
    virtual void WriteH5Params(HighFive::File& file);
    
    //! Write coefficient data in H5
    virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count);
    
  public:
    
    //! Constructor
    SphCoefs(bool verbose=false) : Coefs("sphere", verbose) {}
    
    //! Native file constructor
    SphCoefs(const std::string& file, int stride=1,
	     double tmin=-std::numeric_limits<double>::max(),
	     double tmax= std::numeric_limits<double>::max(),
	     bool verbose=false) :
      Coefs("sphere", verbose)
    {
      readNativeCoefs(file, stride, tmin, tmax);
    }

    //! Copy constructor
    SphCoefs(SphCoefs& p) : Coefs(p)
    {
      coefs = p.coefs;
    }

    //! H5 constructor
    SphCoefs(HighFive::File& file, int stride=0,
	     double tmin=-std::numeric_limits<double>::max(),
	     double tmax= std::numeric_limits<double>::max(),
	     bool verbose=false);
    
    //! Clear coefficient container
    virtual void clear() { coefs.clear(); }

    //! Add a coefficient structure to the container
    virtual void add(CoefStrPtr coef);

    //! Get coefficient matrix at given time
    virtual Eigen::MatrixXcd& operator()(double time);
    
    //! Set coefficient matrix at given time
    virtual void setMatrix(double time, const Eigen::MatrixXcd& mat);

    //! Get coefficient structure at a given time
    virtual std::shared_ptr<CoefStruct> getCoefStruct(double time)
    { return coefs[time]; }

    //! Dump to ascii list for testing
    virtual void dump(int lmin, int lmax, int nmin, int nmax);
    
    //! Get list of coefficient times
    virtual std::vector<double> Times()
    {
      times.clear();
      for (auto t : coefs) times.push_back(t.first);
      return times;
    }
    
    //! Get all coefficients indexed in l, m, n, time
    Eigen::Tensor<std::complex<double>, 3> getAllCoefs();

    /** Get power for the coefficient DB as a function of harmonic
	index.  Time as rows, harmonics as columns.

	@param min is the minimum radial order
	@param max is the maximum radial order
    */
    virtual Eigen::MatrixXd& Power
    (int min=0, int max=std::numeric_limits<int>::max());
    
    //! Make keys for the remaining indices in a subspace
    virtual std::vector<Key> makeKeys(Key k);

    //! Compare two collections of stanzas (this is for testing only)
    virtual bool CompareStanzas(std::shared_ptr<Coefs> check);

    //! Copy all of the data; make new instances of shared pointer
    //! objects
    virtual std::shared_ptr<Coefs> deepcopy();

    //! Zero the existing data
    virtual void zerodata() {
      for (auto v : coefs) v.second->zerodata();
    }

    //! Get lmax
    double lmax() const { return Lmax; }

    //! Get lmax
    double nmax() const { return Nmax; }

  };
  
  /** Derived class for cylindrical coefficients */
  class CylCoefs : public Coefs
  {
  protected:
    
    //! Working array for power
    Eigen::MatrixXd powerE, powerO;
    
    //! Parameters
    int Mmax, Nmax;

    //! The coefficient DB
    std::map<double, CylStrPtr> coefs;

    //! Print module/angle (for diagnostic output)
    bool angle;
    
    //! Read the coefficients
    virtual void readNativeCoefs(const std::string& file,
				 int stride, double tmin, double tmax);

    
    //! Get the YAML config
    virtual std::string getYAML();
    
    //! Write parameter attributes
    virtual void WriteH5Params(HighFive::File& file);
    
    //! Write coefficient data in H5
    virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count);
    
  public:
    
    //! Constructor
    CylCoefs(bool verbose=false) :
      Coefs("cylinder", verbose), angle(false) {}
    
    //! Native file constructor
    CylCoefs(const std::string& file, int stride=1,
	     double tmin=-std::numeric_limits<double>::max(),
	     double tmax= std::numeric_limits<double>::max(),
	     bool verbose=false) :
      Coefs("cylinder", verbose), angle(false)
    {
      readNativeCoefs(file, stride, tmin, tmax);
    }
    
    //! H5 constructor
    CylCoefs(HighFive::File& file, int stride=1,
	     double tmin=-std::numeric_limits<double>::max(),
	     double tmax= std::numeric_limits<double>::max(),
	     bool verbose=false);
    
    //! Copy constructor
    CylCoefs(CylCoefs& p) : Coefs(p)
    {
      coefs = p.coefs;
      angle = p.angle;
    }

    //! Clear coefficient container
    virtual void clear() { coefs.clear(); }

    //! Add a coefficient structure to the container
    virtual void add(CoefStrPtr coef);

    //! Get coefficient matrix at given time
    virtual Eigen::MatrixXcd& operator()(double time);
    
    //! Set coefficient matrix at given time
    virtual void setMatrix(double time, const Eigen::MatrixXcd& mat);
    
    //! Get coefficient structure at a given time
    virtual std::shared_ptr<CoefStruct> getCoefStruct(double time)
    { return coefs[time]; }


    //! Dump to ascii list for testing
    virtual void dump(int mmin, int mmax, int nmin, int nmax);
    
    //! Get list of coefficient times
    virtual std::vector<double> Times()
    {
      times.clear();
      for (auto t : coefs) times.push_back(t.first);
      return times;
    }
    
    //! Get all coefficients indexed in m, n, time
    Eigen::Tensor<std::complex<double>, 3> getAllCoefs();

    /** Get power for the coefficient DB as a function of harmonic
	index.  Time as rows, harmonics as columns.

	@param min is the minimum radial order
	@param max is the maximum radial order
    */
    virtual Eigen::MatrixXd& Power
    (int min=0, int max=std::numeric_limits<int>::max());

    /** Get power for the coefficient DB as a function of harmonic
	index, separating the vertically even and odd symmetric
	coefficients.  Time as rows, harmonics as columns.

	@param min is the minimum radial order
	@param max is the maximum radial order
	@param nodd is the index position of the first vertically odd coefficient
    */
    std::tuple<Eigen::MatrixXd&, Eigen::MatrixXd&> EvenOddPower
    (int nodd=-1, int min=0, int max=std::numeric_limits<int>::max());

    //! Make keys for the remaining indices in a subspace
    virtual std::vector<Key> makeKeys(Key k);

    //! Compare two collections of stanzas (this is for testing only)
    virtual bool CompareStanzas(std::shared_ptr<Coefs> check);

    //! Copy all of the data; make new instances of shared pointer
    //! objects
    virtual std::shared_ptr<Coefs> deepcopy();

    virtual void zerodata() {
      for (auto v : coefs) v.second->zerodata();
    }

    //! Get mmax
    double mmax() const { return Mmax; }

    //! Get lmax
    double nmax() const { return Nmax; }
  };
  

  /** Derived class for a table of data */
  class TableData : public Coefs
  {
  protected:
    
    //! The table DB
    std::map<double, TblStrPtr> coefs;

    //! An alternate packing
    std::vector<std::vector<double>> data;

    //! Read the coefficients
    virtual void readNativeCoefs(const std::string& file,
				 int stride, double tmin, double tmax);
    
    //! Get the YAML config
    virtual std::string getYAML();
    
    //! Write parameter attributes
    virtual void WriteH5Params(HighFive::File& file);
    
    //! Write coefficient data in H5
    virtual unsigned WriteH5Times(HighFive::Group& group, unsigned count);
    
  public:
    
    //! Constructor
    TableData(bool verbose=false) : Coefs("table", verbose) {}

    //! Constructor from vectors
    TableData(const std::vector<double>& time,
	      const std::vector<std::vector<double>>& data,
	      bool verbose=false);

    //! Constructor from file
    TableData(std::string& file, bool verbose=false);
    
    //! H5 constructor
    TableData(HighFive::File& file, int stride=1,
	     double tmin=-std::numeric_limits<double>::max(),
	     double tmax= std::numeric_limits<double>::max(),
	     bool verbose=false);
    
    //! Copy constructor
    TableData(TableData& p) : Coefs(p)
    {
      coefs = p.coefs;
      data  = p.data;
      times = p.times;
    }

    //! Clear coefficient container
    virtual void clear() { coefs.clear(); }

    //! Add a coefficient structure to the container
    virtual void add(CoefStrPtr coef);

    //! Get coefficient matrix at given time
    virtual Eigen::MatrixXcd& operator()(double time);
    
    /** Set coefficient matrix at given time.  This is for pybind11,
	since the operator() will not allow lvalue assignment, it
	seems. */
    virtual void setMatrix(double time, const Eigen::MatrixXcd& mat);

    //! Get coefficient structure at a given time
    virtual std::shared_ptr<CoefStruct> getCoefStruct(double time)
    { return coefs[time]; }

    //! Dump to ascii list for testing
    virtual void dump(int lmin, int lmax, int nmin, int nmax) {}
    
    //! Get list of coefficient times
    virtual std::vector<double> Times() { return times; }
    
    //! Get all coefficients indexed in column, time
    Eigen::MatrixXd getAllCoefs();

    /** Get power for the series.  Not implemented meaningfully
	for table data. */
    virtual Eigen::MatrixXd& Power
    (int min=0, int max=std::numeric_limits<int>::max())
    {
      power.resize(0, 0);
      return power;
    }

    //! Make keys for the remaining indices in a subspace
    virtual std::vector<Key> makeKeys(Key k)
    { return std::vector<Key>(); }

    //! Compare two collections of stanzas (this is for testing only)
    virtual bool CompareStanzas(std::shared_ptr<Coefs> check);

    //! Copy all of the data; make new instances of shared pointer
    //! objects.  This will recursively call the deepcopy() functions
    //! or all dependent objects
    virtual std::shared_ptr<Coefs> deepcopy();

    virtual void zerodata() {
      for (auto v : coefs) v.second->zerodata();
      for (auto & v : data) std::fill(v.begin(), v.end(), 0.0);
    }

  };
  

  using CoefsPtr = std::shared_ptr<Coefs>;
}
// END namespace CoefClasses
  
#endif
  
