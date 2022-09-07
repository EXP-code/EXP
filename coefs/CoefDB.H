#ifndef _CoefDB_H
#define _CoefDB_H

#include <memory>
#include <tuple>
#include <map>
#include <set>

#include <Eigen/Dense>
#include <CoefContainer.H>

namespace MSSA
{
  //! Coefficient and keys
  using CoefElements = std::map<std::string, std::vector<std::vector<int>>>;

  //! An index key
  using Key = std::vector<unsigned>;
  
  //! Coefficient and channel config map/dictionary
  using mssaConfig =
    std::map<std::string, std::tuple<Coefs::CoefsPtr, std::vector<Key>, std::vector<Key>>>;

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

    //! Complex keys
    bool complexKey;

    //! Coefficients
    Coefs::CoefsPtr coefs, coefs1;

    //! Packing routines
    //@{
    //! Reflect type
    void pack_channels();

    //@{
    //! Spherical coefficients
    void pack_sphere();
    void unpack_sphere();

    //@{
    //! Cylindrical coefficients
    void pack_cylinder();
    void unpack_cylinder();
    //@}

    //@{
    //! Table coefficients
    void pack_table();
    void unpack_table();
    //@}
    //@}

    //@{
    //! Not sure we need this since data stream assembly will happen
    //! upstream
    int stride;
    double tmin, tmax;
    //@}

  public:
    //! Harmonics
    using LMKey = std::tuple<unsigned, unsigned>;

    //! Database itself
    std::map<Key, std::vector<double>> data;

    //! Database update
    std::map<Key, std::vector<double>> data1;
    
    //! Times for each slice
    std::vector<double> times;
    
    //! Constructor
    CoefDB(const std::string& name, Coefs::CoefsPtr coefs,
	   const std::vector<Key>& keys, const std::vector<Key>& bkeys, 
	   unsigned index, unsigned stride=1,
	   double tmin=-std::numeric_limits<double>::max(),
	   double tmax= std::numeric_limits<double>::max());


    //! Update working copyof coefficients
    void update(bool zero);

    //! Update and return working copy of coefficients
    Coefs::CoefsPtr endUpdate();
  };

  /**
     A uniform interface to a group of CoefDB instances

     The container reads a YAML configuration file and creates a
     container of coefficient files, one for each stanza in the
     'components' sequence in the YAML file.  There are 4 global
     parameters as well:

     @tmin:   the minimum time to include in the time series
     @tmax:   the maximum time to include in the time series
     @stride: the number of intermediate coefficient bunches to skip
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

    std::vector<Key> keylist;
    std::vector<CoefDBptr> comps;
    std::map<std::string, int> namemap;

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

    //! Null constructor
    CoefContainer() {}

    //! Instantiate all components DBs from YAML spec
    CoefContainer(const std::string spec);
    
    //! Instantiate all components DBs from structure
    CoefContainer(const mssaConfig& config, const std::string spec="");
    
    //! Copy constructor
    CoefContainer(const CoefContainer& p);

    //! Get the pooled key list
    const std::vector<Key>& getKeys() { return keylist; }

    //! Return a particular data stream
    const std::vector<double>& getData(std::vector<unsigned> key)
    {
      unsigned c = key.back();
      key.pop_back();
      if (comps[c]->data.find(key) == comps[c]->data.end()) {
	std::cout << "Keys in comp[" << c << "] :";
	for (auto v : comps[c]->data) {
	  std::cout << " (";
	  for (auto l : v.first) std::cout << l << "|";
	  std::cout << ")";
	}
	std::cout << std::endl;
	std::cout << "Desired key: ";
	std::cout << " (";
	for (auto l : key) std::cout << l << "|";
	std::cout << ")" << std::endl;
	throw std::runtime_error("CoefDB::getData: desired key is not data");
      }
      return comps[c]->data[key];
    }

    //! Copy original data coefficients to working version
    void beginUpdate(bool zero)
    {
      for (auto & v : comps) {
	v->data1 = v->data;
	if (zero) {
	  for (auto & u : v->data1) {
	    std::fill(u.second.begin(), u.second.end(), 0.0);
	  }
	}
      }
    }

    //! Update data in a particular stream
    void setData(std::vector<unsigned> key, const std::vector<double>& data)
    {
      unsigned c = key.back();
      key.pop_back();
      comps[c]->data1[key] = data;
    }

    //! Update the working data to the working coefficients and return a map
    std::map<std::string, Coefs::CoefsPtr> endUpdate()
    {
      std::map<std::string, Coefs::CoefsPtr> ret;
      for (auto & v : comps) {
	ret[v->name] = v->endUpdate();
      }
      return ret;
    }

    //! Are keys complex?  Used to reflect the index type.
    bool isComplex(const std::string& name)
    {
      auto it = namemap.find(name);
      if (it == namemap.end()) {
	throw std::runtime_error
	  ("CoefDB::isComplex: cannot find component with given name");
      }

      return comps[it->second]->complexKey;
    }

    //! Get index with associated name
    int index(const std::string name)
    {
      auto it = namemap.find(name);
      if (it == namemap.end()) return -1;
      else                     return it->second;
    }

    //! Get available names in this DB
    std::vector<std::string> getNames()
    {
      std::vector<std::string> ret;
      for (auto v : namemap) ret.push_back(v.first);
      return ret;
    }

  };

}
// END namespace MSSA

#endif

