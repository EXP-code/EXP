#ifndef _NTC_H
#define _NTC_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <deque>
#include <map>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>

#include <boost/filesystem.hpp>

#include <boost/version.hpp>

#include "pHOT_types.H"
#include "Species.H"
#include "Quantile.H"

// Serialize std::array and std::tuple(...)
//
namespace boost {
    
  namespace serialization {
      
#if BOOST_VERSION <= 105600
    template<class Archive, class T, size_t N>
    void serialize(Archive & ar, std::array<T,N> & a, 
		   const unsigned int version)
    {
      ar & boost::serialization::make_array(a.data(), a.size());
    }
#endif
  
    template<typename Archive, typename T, typename U, typename V>
    void serialize(Archive & ar, std::tuple<T, U, V> & t, 
		   const unsigned int version) 
    {
      ar & boost::serialization::make_nvp("tup0", std::get<0>(t));
      ar & boost::serialization::make_nvp("tup1", std::get<1>(t));
      ar & boost::serialization::make_nvp("tup2", std::get<2>(t));
    }

  } // namespace serialization

} // namespace boost


template<typename T, typename U, typename V>
std::ostream& operator<< (std::ostream& out, const std::tuple<T, U, V>& t);

std::ostream& operator<< (std::ostream& out, const sKeyPair& t);


namespace NTC {

  //! Container for interaction list to add UI to data types
  class Interact
  {

  public:

    //! Type of particle in interaction
    enum pType {simple, neutral, ion, electron};

    //! Interaction element: interaction species and charge
    typedef std::pair<pType, unsigned short> pElem;

    //! Interaction type with interaction elements
    typedef std::tuple<unsigned short, pElem, pElem> T;

    //! Default null pElem
    static pElem const pdef;

    //! Default electron pElem
    static pElem const edef;

    //! Default interaction key for singleton
    static T single;

  private:

    //! For labeling output
    static std::map<pType, std::string> pLabs;

  public:

    //! Array-like data implemented by a map
    std::map<T, double> v;

    //! Single double value
    double d;

    //! Constructor (for initialization)
    Interact() : d(0) {}

    //! Is map empty?
    bool operator!() { return v.size() == 0; }

    //! Get double from the map
    double& operator[](T t) { return v[t]; }

    //! Get single double
    double& operator()() { return d; }

    //! Return interaction label
    static const std::string& label(const pElem& t) { return pLabs[t.first]; }

  };

  class NTCitem
  {
    friend class boost::serialization::access;
    friend class NTCdb;
    
  public:

    typedef NTC::Interact::T T;

    //! Default interaction key for singleton
    static T single;

    //@{
    //! NTC structures

    class udMap : public std::map<T, double>
    {
      friend class boost::serialization::access;

    public:

      udMap() 
      {
	(*this)[single] = 0.0;
      }

      udMap(T t, double z)
      {
	(*this)[t] = z;
      }

    private:

      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
      {
	typedef std::map<T, double> Tup3DoubleMap;
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Tup3DoubleMap);
      }
      
    };

    class uqMap : public std::map<T, Quantile>
    {
      friend class boost::serialization::access;

    public:

      uqMap() {}

      uqMap(const Quantile& q) 
      {
	(*this)[single] = Quantile();
	(*this)[single] = q;
      }

      uqMap(T t, const Quantile& q)
      {
	(*this)[t] = Quantile();
	(*this)[t] = q;
      }

    private:

      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
      {
	typedef std::map<T, Quantile> Tup3QuantileMap;
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Tup3QuantileMap);
      }
    };

    typedef std::map<sKeyPair, udMap> vcMap;
    typedef std::map<sKeyPair, uqMap> qpMap;
    //@}
    
  protected:
    
    //@{
    //! Quantile init data
    
    //@{
    //! The desired quantiles.  You can add as many of these 
    //! as you want.
    const int Nequal = 10;
    const double u1 = 0.01, u2 = 0.05, u3 = 0.95, u4 = 0.99;
    const std::vector<double> qs = {u1, u2, u3, u4};
    //@}
    
    //! Main NTC data structures
    qpMap db;
    
    void debug();

    //@{
    //! For debugging
    static unsigned maxrpt;
    static unsigned skprpt;
    key_type caller;
    unsigned report;
    bool inTest;
    
    void Test();
    void setKey(const key_type& k) { caller = k; }
    //@}
    
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(db);
      ar & BOOST_SERIALIZATION_NVP(caller);
      ar & BOOST_SERIALIZATION_NVP(report);
      ar & BOOST_SERIALIZATION_NVP(inTest);
    }
    
  private:

    static unsigned instance;

    void oddBall(sKeyPair indx, const T& intr)
    {
      speciesKey k1 = indx.first, k2 = indx.second;
      if (
	  (k1.first != 0 or k1.second != 0 or
	   k2.first != 0 or k2.second != 0) and intr == single) 
	{
	  std::cout << "Odd ball: " << indx
		    << " : " << std::get<0>(intr) << std::endl;
	}
    }

  public:
    
    //! Constructor
    NTCitem() 
    { 
      report = 0; 
      inTest = false;

      instance++;
    }
    
    //! Copy constructor
    NTCitem(const NTCitem& p) 
    {
      db.insert(p.db.begin(), p.db.end());
      report = 0;
      inTest = false;

      instance++;
    }
    
    //! Assignment operator
    NTCitem &operator=(const NTCitem& p)
    {
      db.insert(p.db.begin(), p.db.end());
      report = 0;
      inTest = false;

      instance++;

      return *this;
    }

    //! Destructor (for debugging live instance count)
    ~NTCitem()
    {
      instance--;
    }
    
    //! Threshold <cross section> x <relative velocity> value
    static double Min;
    
    //! Default <cross section> x <relative velocity> value
    static double Def;
    
    //@{
    //! Return quantile value for a single index
    double Prob(sKeyPair indx, const T& intr, double x);
    double Prob(sKeyPair indx, double x) 
    {
      return Prob(indx, single, x);
    }
    //@}
    
    //! Return quantile values for all indices
    vcMap Prob(double x);
    
    //@{
    //! Return value for quantile q for a single index
    double CrsVel(sKeyPair indx, double p) 
    {
      return CrsVel(indx, single, p);
    }
    double CrsVel(sKeyPair indx, const T& intr, double p);
    //@}
    
    //! Return values for quantile q for all indices
    vcMap CrsVel(double p);
    
    //@{
    //! Add a single value to the list
    void  Add(sKeyPair indx, const T& intr, double x);
    void  Add(sKeyPair indx, double x) { Add(indx, single, x); }
    //@}
    
    //@{
    //! Ready for estimation owing to a full set of points
    bool Ready(sKeyPair indx, const T& intr)
    {
      oddBall(indx, intr);

      qpMap::iterator it = db.find(indx);
      if (it == db.end()) return false;
      else {
	uqMap::iterator jt = it->second.find(intr);
	if (jt == it->second.end()) return false;
	return jt->second.full();
      }
    }
    bool Ready(sKeyPair indx) { return Ready(indx, single); }
    //@}

    //! Node sends to root
    void send()
    {
      unsigned long clr = caller;
      unsigned siz      = db.size();
      
      MPI_Send(&clr, 1, MPI_UNSIGNED_LONG, 0, 1110, MPI_COMM_WORLD);
      MPI_Send(&siz, 1, MPI_UNSIGNED,      0, 1111, MPI_COMM_WORLD);

      for (auto &p : db) {
	speciesKey p1 = p.first.first;
	speciesKey p2 = p.first.second;
	
	MPI_Send(&p1.first,  1, MPI_UNSIGNED_SHORT, 0, 1112, MPI_COMM_WORLD);
	MPI_Send(&p1.second, 1, MPI_UNSIGNED_SHORT, 0, 1113, MPI_COMM_WORLD);
	MPI_Send(&p2.first,  1, MPI_UNSIGNED_SHORT, 0, 1114, MPI_COMM_WORLD);
	MPI_Send(&p2.second, 1, MPI_UNSIGNED_SHORT, 0, 1115, MPI_COMM_WORLD);
	
	unsigned siz2 = p.second.size();
	MPI_Send(&siz2, 1, MPI_UNSIGNED,   0, 1116, MPI_COMM_WORLD);

	unsigned short u;	// temporary

	for (auto &q : p.second) {
	  MPI_Send(&(u=std::get<0>(q.first)), 1, MPI_UNSIGNED_SHORT, 0, 1117, MPI_COMM_WORLD);
	  MPI_Send(&(u=std::get<1>(q.first).first),  1, MPI_UNSIGNED_SHORT, 0, 1118, MPI_COMM_WORLD);
	  MPI_Send(&(u=std::get<1>(q.first).second), 1, MPI_UNSIGNED_SHORT, 0, 1119, MPI_COMM_WORLD);
	  MPI_Send(&(u=std::get<2>(q.first).first),  1, MPI_UNSIGNED_SHORT, 0, 1120, MPI_COMM_WORLD);
	  MPI_Send(&(u=std::get<2>(q.first).second), 1, MPI_UNSIGNED_SHORT, 0, 1121, MPI_COMM_WORLD);

	  q.second.send();
	}
      }
      
    }
    
    // Root node initializes from node
    //
    void recv(int id)
    {
      MPI_Recv(&caller, 1, MPI_UNSIGNED_LONG, id, 1110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      unsigned siz;
      
      MPI_Recv(&siz, 1, MPI_UNSIGNED, id, 1111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      for (unsigned k=0; k<siz; k++) {
	
	speciesKey p1, p2;
	
	MPI_Recv(&p1.first,  1, MPI_UNSIGNED_SHORT, id, 1112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&p1.second, 1, MPI_UNSIGNED_SHORT, id, 1113, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	MPI_Recv(&p2.first,  1, MPI_UNSIGNED_SHORT, id, 1114, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&p2.second, 1, MPI_UNSIGNED_SHORT, id, 1115, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	unsigned siz2;
      
	MPI_Recv(&siz2, 1, MPI_UNSIGNED, id, 1116, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	sKeyPair key(p1, p2);

	for (unsigned j=0; j<siz2; j++) {
	  T t;
	  MPI_Recv(&(get<0>(t)),        1, MPI_UNSIGNED_SHORT, id, 1117, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  MPI_Recv(&(get<1>(t).first),  1, MPI_UNSIGNED_SHORT, id, 1118, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&(get<1>(t).second), 1, MPI_UNSIGNED_SHORT, id, 1119, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  MPI_Recv(&(get<2>(t).first),  1, MPI_UNSIGNED_SHORT, id, 1120, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&(get<2>(t).second), 1, MPI_UNSIGNED_SHORT, id, 1121, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  db[key][t].recv(id);

	} // db2 size loop

      } // db1 size loop

} // function end

    static unsigned live() { return instance; }
    
  };
  
  typedef std::map<key_type, NTCitem> NTCdata;
  
  
  //! Provide a nice interface for the NTCdata std::map
  class NTCdb
  {
    friend class boost::serialization::access;
    
  protected:
    
    //! The std::map
    NTCdata data;
    
    //! Current serialization time
    double curr_time;

    //! Next serialization time
    double next_time;
    
    //! Roundoff offset
    const double eps = 1.0e-8;

    //! Persistent copy
    NTCdata pdata;
    
    //! Get copies from all nodes to save
    void sync() 
    {
      // Insert root's nodes into pdata
      if (myid==0) pdata.insert(data.begin(), data.end());
      
      unsigned siz;
#ifdef INT128
      uint128 k;
      unsigned long k1, k2;
#else
      unsigned long k;
#endif
      
      for (int i=1; i<numprocs; i++) {
	if (i==myid) {
	  siz = data.size();
	  MPI_Send(&siz, 1, MPI_UNSIGNED, 0, 1110, MPI_COMM_WORLD);
	  for (auto &p : data) {
	    k = p.first;
#ifdef INT128
	    k.decompose(k1, k2);
	    MPI_Send(&k1, 1, MPI_UNSIGNED_LONG, 0, 1111, MPI_COMM_WORLD);
	    MPI_Send(&k2, 1, MPI_UNSIGNED_LONG, 0, 1112, MPI_COMM_WORLD);
#else
	    MPI_Send(&k,  1, MPI_UNSIGNED_LONG, 0, 1111, MPI_COMM_WORLD);
#endif
	    p.second.send();
	  }
	  
	}
	
	if (0==myid) {
	  MPI_Recv(&siz, 1, MPI_UNSIGNED, i, 1110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for (unsigned j=0; j<siz; j++) {
#ifdef INT128
	    MPI_Recv(&k1, 1, MPI_UNSIGNED_LONG, i, 1111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(&k2, 1, MPI_UNSIGNED_LONG, i, 1112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    k.recompose(k1, k2);
#else
	    MPI_Recv(&k,  1, MPI_UNSIGNED_LONG, i, 1111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
	    
	    // Get the instance from the other node
	    //
	    pdata[k].recv(i);
	  }
	}
	// Need this barrier to prevent the MPI buffer cache from
	// getting insanely large
	//
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    
    //! Save contents to a file
    void save_myself(const std::string& filename)
    {
      // Synchronize
      sync();
      
      // Only the root process makes an archive
      if (myid==0) {
	std::ofstream ofs;
	bool exists = boost::filesystem::exists(filename);
	std::string tfile(".tmp." + filename);
	
	if (exists) {
	  ofs.open(tfile.c_str());
	} else {
	  ofs.open(filename.c_str());
	}

	boost::archive::xml_oarchive oa(ofs);
	oa << BOOST_SERIALIZATION_NVP(pdata);
	oa << BOOST_SERIALIZATION_NVP(curr_time);
	oa << BOOST_SERIALIZATION_NVP(next_time);

	if (exists) {
	  boost::filesystem::rename(tfile, filename);
	}

	if (chatty) {
	  std::cout << std::endl << std::string(70, '-') << std::endl
		    << "-- Took NTC checkpoint at time <" << curr_time << ">" << std::endl
		    << "-- Next NTC checkpoint at time <" << next_time << ">" << std::endl
		    << std::string(70, '-') << std::endl;
	}
      }
      
      // Clean up node-aggregated db
      pdata.clear();

      // Debugging
      if (myid==0) {
	std::cout << std::string(40, '-') << std::endl
		  << "--- Quantile instances: " << Quantile::live()
		  << std::endl
		  << "--- NTCitem instances:  " << NTCitem::live()
		  << std::endl
		  << std::string(40, '-') << std::endl;
      }
    }
    
    //! Restore contents from a file
    void restore_myself(const std::string& filename)
    {
      // open the archive; all processes
      std::ifstream ifs(filename.c_str());
      boost::archive::xml_iarchive ia(ifs);
      
      // restore the data from the archive
      ia >> BOOST_SERIALIZATION_NVP(data);
      ia >> BOOST_SERIALIZATION_NVP(curr_time);
      ia >> BOOST_SERIALIZATION_NVP(next_time);

      if (myid==0 and chatty) {
	std::cout << std::endl << std::string(70, '-') << std::endl
		  << "-- Restoring NTC checkpoint from time <" << curr_time << ">" << std::endl
		  << "-- Current simulation time is <" << tnow << ">" << std::endl
		  << "-- Next NTC checkpoint scheduled for <" << next_time << ">" << std::endl
		  << std::string(70, '-') << std::endl;
      }

    }
    
  public:
    
    //! Chatty output for debugging
    static bool chatty;

    //! Save interval
    static unsigned intvl;
    
    //! Constructor
    NTCdb()
    {
      // Attempt to restore from save set
      std::string filename = ".ntcdb." + runtag;
      if ( boost::filesystem::exists(filename) ) {
	restore_myself(filename);
      }
      //
      curr_time = tnow;
      next_time = curr_time + dtime*(static_cast<double>(intvl)-eps);
    }
    
    //! Return NTCitem instance for cell key k
    NTCitem& operator[](const key_type& k);
    
    //! Update the serialization file
    void update()
    {
      if (tnow >= next_time) {
	curr_time = tnow;
	next_time = curr_time + dtime*(static_cast<double>(intvl)-eps);
	std::string filename = ".ntcdb." + runtag;
	save_myself(filename);
      }
    }
    
    //! Serialize to finish
    void finish()
    {
      std::string filename = ".ntcdb." + runtag;
      save_myself(filename);
    }

  };

} // End: namespace NTC

#endif
  
