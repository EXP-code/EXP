#ifndef _CoefContainer_H
#define _CoefContainer_H

#include <memory>
#include <tuple>
#include <map>
#include <set>

#include <Eigen/Dense>
#include <Coefficients.H>

namespace MSSA
{
  //! Coefficient and keys
  using CoefElements = std::map<std::string, std::vector<std::vector<int>>>;

  //! An index key
  using Key = std::vector<unsigned>;
  
  //! Key comparison with length ordering
  struct mSSAkeyCompare
  {
    bool operator()(const std::vector<unsigned>& a,
		    const std::vector<unsigned>& b) const
    {
				// Get component indices
      unsigned A = a.back(), B = b.back();

      if (A < B)		// Order by component
	return true;

      if (A == B) {		// Order by index
	size_t sz = a.size()-1;
	for (size_t i=0; i<sz; i++) {
	  if (a[i] < b[i]) return true;
	  if (a[i] > b[i]) return false;
	}
      }
      return false;
    }
  };

  //! Coefficient and channel config map/dictionary
  using mssaConfig =
    std::map<std::string, std::tuple<CoefClasses::CoefsPtr, std::vector<Key>, std::vector<Key>>>;

  //! Detrending style
  enum class TrendType {perChannel, totVar, totPow};

  //! Get enum names from a constant map (this is tedious; there is
  //! probably an easier way)
  const std::map<int, std::string> getTrendType
    { {static_cast<std::underlying_type<TrendType>::type>(TrendType::perChannel), "per channel"   },
      {static_cast<std::underlying_type<TrendType>::type>(TrendType::totVar),     "total variance"},
      {static_cast<std::underlying_type<TrendType>::type>(TrendType::totPow),     "total power"   } };


  class CoefContainer;

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
    std::vector<Key> keys, keys0;

    //! List of keys for background
    std::vector<Key> bkeys, bkeys0;

    //! Complex keys
    bool complexKey;

    //! Coefficients
    CoefClasses::CoefsPtr coefs;

    //! Packing routines
    //@{
    //! Reflect type
    void pack_channels();

    //@{
    //! Spherical coefficients
    void pack_sphere();
    void unpack_sphere();
    void restore_background_sphere();
    //@}

    //@{
    //! Cylindrical coefficients
    void pack_cylinder();
    void unpack_cylinder();
    void restore_background_cylinder();
    //@}

    //@{
    //! Spherical field coefficients
    void pack_sphfld();
    void unpack_sphfld();
    //@}

    //@{
    //! Cylindrical coefficients
    void pack_cylfld();
    void unpack_cylfld();
    //@}

    //@{
    //! Slab coefficients
    void pack_slab();
    void unpack_slab();
    //@}

    //@{
    //! Cube coefficients
    void pack_cube();
    void unpack_cube();
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
    std::map<Key, std::vector<double>, mSSAkeyCompare> data;

    //! Times for each slice
    std::vector<double> times;
    
    //! Null constructor
    CoefDB() {}

    //! Constructor
    CoefDB(const std::string& name, CoefClasses::CoefsPtr coefs,
	   const std::vector<Key>& keys, const std::vector<Key>& bkeys, 
	   unsigned index, unsigned stride=1,
	   double tmin=-std::numeric_limits<double>::max(),
	   double tmax= std::numeric_limits<double>::max());

    //! Update and return working copy of coefficients
    CoefClasses::CoefsPtr endUpdate();

    //! Copy all of the data; make new instances of shared pointer
    //! objects
    std::shared_ptr<CoefDB> deepcopy();

    //! Copy data from backround back to the coefficient database
    void background();

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

    //! Copy all of the data; make new instances of shared pointer
    //! objects
    std::shared_ptr<CoefContainer> deepcopy();

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

    //! Update data in a particular stream
    void setData(std::vector<unsigned> key, const std::vector<double>& data)
    {
      unsigned c = key.back();
      key.pop_back();
      comps[c]->data[key] = data;
    }

    //! Copy data from backround back to the coefficient database
    void background()
    {
      for (auto v : comps) v->background();
    }

    //! Update the working data to the working coefficients and return a map
    std::map<std::string, CoefClasses::CoefsPtr> endUpdate()
    {
      std::map<std::string, CoefClasses::CoefsPtr> ret;
      for (auto & v : comps) {
	ret[v->name] = v->endUpdate();
      }
      return ret;
    }

    //! Are keys complex?  Used by clients to reflect the index type.
    bool isComplex(const std::string& name)
    {
      auto it = namemap.find(name);
      if (it == namemap.end()) {
	throw std::runtime_error
	  ("CoefDB::isComplex: cannot find component with given name");
      }

      return comps[it->second]->complexKey;
    }

    //! Get index with associated name.  Returns -1 if name is not found.
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
  // END CoefContainer

}
// END namespace MSSA

#endif

