#ifndef _CoefDB_H
#define _CoefDB_H

#include <memory>
#include <Eigen/Dense>
#include "Coefs.H"

namespace CoefDB {

  //! An index key
  using Key = std::vector<unsigned>;
  
  //! Detrending style
  enum class TrendType {perChannel, totVar, totPow};

  //! Base class for per component coefficients
  class CoefDB
  {
    friend class CoefContainer;

  protected:
    //! Name of component
    std::string name;

    //! Index of component
    unsigned index;

    //! List of desired keys
    std::vector<Key> keys;

    //! List of keys for background
    std::vector<Key> bkeys;

  public:
    //! Harmonics
    using LMKey = std::tuple<unsigned, unsigned>;

    //! Database itself
    std::map<Key, std::vector<double>> coefs;
    
    //! Times for each slice
    std::vector<double> times;
    
    //! Constructor
    CoefDB(const std::string& name, YAML::Node& node, unsigned index);

    //! Read coefficients from EXP file
    virtual void read(const std::string& file, unsigned stride, double tmin, double tmax) = 0;

    //! Write reconstructions for EXP rendering
    virtual void
    write(int ncomp, const TrendType& type, bool useMean,
	  double totPow, double totVar,
	  std::map<Key, double>& mean,
	  std::map<Key, double>& var,
	  std::map<Key, Eigen::MatrixXd>& RC,
	  const std::vector<std::set<int>>& groups,
	  const std::vector<double>& times,
	  const std::string& runtag) = 0;

    //! Write reconstructions in ascii
    virtual void
    write_ascii(int ncomp, const TrendType& type, bool useMean,
		double totPow, double totVar,
		std::map<Key, double>& mean,
		std::map<Key, double>& var,
		std::map<Key, Eigen::MatrixXd>& RC,
		const std::vector<std::set<int>>& groups,
		const std::vector<double>& times,
		const std::string& runtag) = 0;
  };

  class CylDB : public CoefDB
  {
  protected:

    //! Cache of original data
    std::map<double, CylCoefsPtr> data;

    //! Read coefficients from EXP file
    void read(const std::string& file, unsigned stride=1,
	      double tmin=-std::numeric_limits<double>::max(),
	      double tmax= std::numeric_limits<double>::max());

    int mmax, nmax, ntimes;

  public:
    
    //! Constructor
    CylDB(const std::string& name, YAML::Node& node, unsigned index, unsigned stride=1,
	  double tmin=-std::numeric_limits<double>::max(),
	  double tmax= std::numeric_limits<double>::max());

    //! Write reconstructions for EXP rendering
    void write(int ncomp, const TrendType& type, bool useMean,
	       double totPow, double totVar,
	       std::map<Key, double>& mean,
	       std::map<Key, double>& var,
	       std::map<Key, Eigen::MatrixXd>& RC,
	       const std::vector<std::set<int>>& groups,
	       const std::vector<double>& times,
	       const std::string& runtag);

    //! Write reconstructions in ascii
    void write_ascii(int ncomp, const TrendType& type, bool useMean,
		     double totPow, double totVar,
		     std::map<Key, double>& mean,
		     std::map<Key, double>& var,
		     std::map<Key, Eigen::MatrixXd>& RC,
		     const std::vector<std::set<int>>& groups,
		     const std::vector<double>& times,
		     const std::string& runtag);

  };

  class SphDB : public CoefDB
  {
  protected:

    //! Cache of original data
    std::map<double, SphCoefsPtr> data;

    //! Read coefficients from EXP file
    void read(const std::string& file, unsigned stride=1,
	      double tmin=-std::numeric_limits<double>::max(),
	      double tmax= std::numeric_limits<double>::max());
      
    int lmax, nmax, ntimes;

  public:
    
    SphDB(const std::string& name, YAML::Node& node, unsigned index, unsigned stride=1,
	  double tmin=-std::numeric_limits<double>::max(),
	  double tmax= std::numeric_limits<double>::max());
  
    //! Write reconstructions for EXP rendering
    void write(int ncomp, const TrendType& type, bool useMean,
	       double totPow, double totVar,
	       std::map<Key, double>& mean,
	       std::map<Key, double>& var,
	       std::map<Key, Eigen::MatrixXd>& RC,
	       const std::vector<std::set<int>>& groups,
	       const std::vector<double>& times,
	       const std::string& runtag);


    //! Write reconstructions in ascii
    void write_ascii(int ncomp, const TrendType& type, bool useMean,
		     double totPow, double totVar,
		     std::map<Key, double>& mean,
		     std::map<Key, double>& var,
		     std::map<Key, Eigen::MatrixXd>& RC,
		     const std::vector<std::set<int>>& groups,
		     const std::vector<double>& times,
		     const std::string& runtag);
  };

  class TableDB : public CoefDB
  {
  protected:

    //! Cache of original data
    std::map<double, std::vector<double>> data;

    //! Read coefficients from EXP file
    void read(const std::string& file, unsigned stride=1,
	      double tmin=-std::numeric_limits<double>::max(),
	      double tmax= std::numeric_limits<double>::max());
      
    int nfields, ntimes;

  public:
    
    TableDB(const std::string& name, YAML::Node& node, unsigned index, unsigned stride=1,
	    double tmin=-std::numeric_limits<double>::max(),
	    double tmax= std::numeric_limits<double>::max());
  
    //! Binary write reconstruction
    void write(int ncomp, const TrendType& type, bool useMean,
	       double totPow, double totVar,
	       std::map<Key, double>& mean,
	       std::map<Key, double>& var,
	       std::map<Key, Eigen::MatrixXd>& RC,
	       const std::vector<std::set<int>>& groups,
	       const std::vector<double>& times,
	       const std::string& runtag)
    {
      // Do nothing
    }

    //! Write reconstructions in ascii
    void write_ascii(int ncomp, const TrendType& type, bool useMean,
		     double totPow, double totVar,
		     std::map<Key, double>& mean,
		     std::map<Key, double>& var,
		     std::map<Key, Eigen::MatrixXd>& RC,
		     const std::vector<std::set<int>>& groups,
		     const std::vector<double>& times,
		     const std::string& runtag);
  };


  /**
     A uniform interface to a group of CoefDB instances

     The container reads a YAML configuration file and creates a
     container of coefficient files, one for each stanza in the
     'components' sequence in the YAML file.  There are 4 global
     parameters as well:

     @tmin:   the minimum time to include in the time series
     @tmax:   the maximum time to include in the time series
     @strid:  the number of intermediate coefficient bunches to skip
     @runtag: an internal tag for prefixing output files

     An example YAML file might look something like this:

---
tmin: 0.0
tmax: 2.0
runtag: run1
components:
  disk:
    filename: outcoef.star.run1
    geometry: cylinder
    keys:
      - [0, 16, 0]
      - [0, 17, 0]
      - [0, 18, 0]
      - [0, 19, 0]
  halo:
    filename: outcoef.dark.run1
    geometry: sphere
    keys:
      - [1, 0, 4, 0]
      - [1, 0, 5, 0]
      - [1, 0, 6, 0]
      - [1, 0, 7, 0]
...

     This reads in two coefficient files, outcoef.star.run1 and
     outcoef.halo.run1, respectively.  The first is a cylindrical disk
     and the second is a spherical halo.  The coefficient type is
     specified by the 'geometry' tag.

     The 'keys' describe the coefficient channels from each set to be
     included in the MSSA analysis.  

     In the cylindrical case, the key index sequence is: azimuthal
     harmonic (in this case: m=0), the radial orders (in this case:
     n=16,17,18,19), and the cosine (0) or sine (1) coefficients.  In
     this case, only cosines.

     In the spherical case, the key index sequence is: the spherical
     harmonic pair l and m (in this case: l=1, m=0), the radial orders
     (in this case: n=4,5,6,7), and the cosine (0) or sine (1)
     coefficients.  In this case, only cosines.
   */

  class CoefContainer
  {
  protected:
    using CoefDBptr = std::shared_ptr<CoefDB>;

    std::vector<CoefDBptr> comps;
    std::vector<Key> keylist;

  public:
    
    //@{
    //! Time limits (-inf, inf, by default)
    double tmin, tmax;
    //@}

    //! Time slice stride (1 by default)
    unsigned stride;

    //! Time array
    std::vector<double> times;

    //! Run tag for output
    std::string runtag;

    //! Instantiate all components DBs from YAML spec
    CoefContainer(const std::string& spec);
    
    //! Get the pooled key list
    const std::vector<Key>& getKeys() { return keylist; }

    //! Return a particular data stream
    const std::vector<double>& getData(std::vector<unsigned> key)
    {
      unsigned c = key.back();
      key.pop_back();
      return comps[c]->coefs[key];
    }

    //! Write all reconstructions
    void write_recon(int ncomp, const TrendType& type, bool useMean,
		     double totPow, double totVar,
		     std::map<Key, double>& mean,
		     std::map<Key, double>& var,
		     std::map<Key, Eigen::MatrixXd>& RC,
		     const std::vector<std::set<int>>& groups)
    {
      for (auto c : comps) {
	c->write(ncomp, type, useMean, totPow, totVar,
		 mean, var, RC, groups, times, runtag);
      }
    }

    //! Write all to ascii 
    void write_ascii(int ncomp, const TrendType& type, bool useMean,
		     double totPow, double totVar,
		     std::map<Key, double>& mean,
		     std::map<Key, double>& var,
		     std::map<Key, Eigen::MatrixXd>& RC,
		     const std::vector<std::set<int>>& groups)
    {
      for (auto c : comps) {
	c->write_ascii(ncomp, type, useMean, totPow, totVar,
		       mean, var, RC, groups, times, runtag);
      }
    }

  };

}

#endif

