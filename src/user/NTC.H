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

#include "pHOT_types.H"
#include "Species.H"
#include "Quantile.H"

//
// Tuple arithmetic
//

//----------------------------------------------------------------------
// Add two tuples
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator+(const std::tuple<T1, T2, T3>& a, 
				 const std::tuple<T1, T2, T3>& b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) + std::get<0>(b), 
			      std::get<1>(a) + std::get<1>(b), 
			      std::get<2>(a) + std::get<2>(b) );
  return ret;
};

//----------------------------------------------------------------------
// Subtract two tuples
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator-(const std::tuple<T1, T2, T3>& a, 
				 const std::tuple<T1, T2, T3>& b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) - std::get<0>(b), 
			      std::get<1>(a) - std::get<1>(b), 
			      std::get<2>(a) - std::get<2>(b) );
  return ret;
};

//----------------------------------------------------------------------
// Multiply a tuple by a scalar
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S>
std::tuple<T1, T2, T3> operator*(const std::tuple<T1, T2, T3>& a, const S b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) * b,
			      std::get<1>(a) * b,
			      std::get<2>(a) * b );
  return ret;
};

//----------------------------------------------------------------------
// Divide a tuple by a scalar
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S>
std::tuple<T1, T2, T3> operator/(const std::tuple<T1, T2, T3>& a, const S b)
{
  std::tuple<T1, T2, T3> ret( std::get<0>(a) / b,
			      std::get<1>(a) / b,
			      std::get<2>(a) / b );
  return ret;
};

//----------------------------------------------------------------------
// Add two tuples (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator+=(std::tuple<T1, T2, T3>& a, const std::tuple<T1, T2, T3>& b)
{
  std::get<0>(a)  += std::get<0>(b);
  std::get<1>(a)  += std::get<1>(b);
  std::get<2>(a)  += std::get<2>(b);

  return a;
};

//----------------------------------------------------------------------
// Multiple a tuple by a scalar (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S >
std::tuple<T1, T2, T3> operator*=(std::tuple<T1, T2, T3>& a, const S b)
{
  std::get<0>(a) *= b;
  std::get<1>(a) *= b;
  std::get<2>(a) *= b;

  return a;
};

//----------------------------------------------------------------------
// Divide a tuple by a scalar (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3, class S >
std::tuple<T1, T2, T3> operator/=(std::tuple<T1, T2, T3>& a, const S b)
{
  std::get<0>(a) /= b;
  std::get<1>(a) /= b;
  std::get<2>(a) /= b;

  return a;
};

//----------------------------------------------------------------------
// Subtract two tuples (non const version)
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::tuple<T1, T2, T3> operator-=(std::tuple<T1, T2, T3>& a, const std::tuple<T1, T2, T3>& b)
{
  std::get<0>(a) -= std::get<0>(b);
  std::get<1>(a) -= std::get<1>(b);
  std::get<2>(a) -= std::get<2>(b);

  return a;
};

//----------------------------------------------------------------------
// iostream manipulators
//----------------------------------------------------------------------

template< class T1, class T2, class T3 >
std::ostream& operator<< (std::ostream& stream, 
			  const std::tuple<T1, T2, T3>& v)
{
  stream << "["  << std::setw(12) << std::get<0>(v) 
	 << ", " << std::setw(12) << std::get<1>(v) 
	 << ", " << std::setw(12) << std::get<2>(v) << "]";

  return stream;
}

//! Syntatic glue class to change output for std::vector<double>
class vbkts 
{
private:
  std::vector<double> data;

public:
  //! Construct the wrapper
  vbkts(std::vector<double>& d) { data = d; }
  
  //! Print to stream
  friend std::ostream& operator<<(std::ostream& stream, const vbkts& v)
  {
    stream << "[";
    for (auto d : v.data) stream << std::setw(12) << d;
    stream << "]";
    return stream;
  }
};

// Use a recursive C++-11 variadic template to serialize std::tuple
//
namespace boost {
    
  namespace serialization {
      
    // The generic recursion function for all but the first (<0>)
    // tuple position
    //
    template<uint N>
    struct Serialize
    {
      template<class Archive, typename... Arguments>
      static void serialize(Archive & ar, std::tuple<Arguments...> & t, 
			    const unsigned int version)
      {
	// Make the id name for the name-value pair
	std::ostringstream s; s << "tuple_" << N-1;
	// Archive
	ar & make_nvp(s.str().c_str(), std::get<N-1>(t));
	// Do next element
	Serialize<N-1>::serialize(ar, t, version);
      }
    };
    
    // Specialize the <0> case to stop the recursion
    //
    template<>
    struct Serialize<0>
    {
      template<class Archive, typename... Arguments>
      static void serialize(Archive & ar, std::tuple<Arguments...> & t, 
			    const unsigned int version)
      {
	// Do nothing, end recursion
      }
    };
    
    // Finally, the serialize template for std::tuple
    template<class Archive, typename... Arguments>
    void serialize(Archive & ar, std::tuple<Arguments...> & t, 
		   const unsigned int version)
    {
      Serialize<sizeof...(Arguments)>::serialize(ar, t, version);
    }
    
  }
}
  

// Serialize std::array
//
namespace boost {
    
  namespace serialization {
      
    template<class Archive, class T, size_t N>
    void serialize(Archive & ar, std::array<T,N> & a, 
		   const unsigned int version)
    {
      ar & boost::serialization::make_array(a.data(), a.size());
    }

  } // namespace serialization
  
} // namespace boost



namespace NTC {

  //! NTC structures defintion
  typedef std::tuple<double, double, double> vcTup;
  
  class QuantileArray
  {
    friend class boost::serialization::access;
    
  public:
    
    typedef std::array<Quantile, 3> Q3;
    
  private:
    std::map<double, Q3> tt;
    
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(tt);
    }
    
  public:
    //! Null constructor
    QuantileArray() {}
    
    //! Construct an array of quantile pointers
    QuantileArray(const std::vector<double>& e)
    {
      initialize(e);
    }
    
    //! Copy constructor
    QuantileArray(const QuantileArray& p)
    {
      tt = p.tt;
    }

    //! Initialize the array
    void initialize(const std::vector<double>& e)
    {
      tt.erase(tt.begin(), tt.end());
      
      for (auto p : e) {
	for (auto &q : tt[p]) q.reset(p);
      }
    }
    
    void add(const vcTup& v) 
    {
      for (auto &t : tt) {
	t.second[0](std::get<0>(v));
	t.second[1](std::get<1>(v));
	t.second[2](std::get<2>(v));
      }
    }
    
    Q3 operator[](double p)
    {
      return tt[p];
    }
    
    //! Node sends to root
    void send()
    {
      // Sends the array of Quantiles
      //
      unsigned siz = tt.size();
      // std::cout << "Node " << myid << " sending " << siz
      // << " quantiles" << std::endl;
      MPI_Send(&siz, 1, MPI_UNSIGNED, 0, 1130, MPI_COMM_WORLD);
      for (auto &p : tt) {
	double d = p.first;
	MPI_Send(&d, 1, MPI_DOUBLE, 0, 1131, MPI_COMM_WORLD);
	p.second[0].send();
	p.second[1].send();
	p.second[2].send();
      }
    }
    
    //! Root initializes itself from node's data
    void recv(int id) 
    {
      unsigned siz;
      MPI_Recv(&siz, 1, MPI_UNSIGNED, id, 1130, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      
      // std::cout << "Root receiving " << siz
      // << " quantiles" << std::endl;

      // Receives the array of quantiles from a node
      //
      for (unsigned j=0; j<siz; j++) {
	double d;
	MPI_Recv(&d, 1, MPI_DOUBLE, id, 1131, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	tt[d][0].recv(id);
	tt[d][1].recv(id);
	tt[d][2].recv(id);
      }
    }
    
  }; // end: QuantileArray
  
  
  class NTCitem
  {
    friend class boost::serialization::access;
    friend class NTCdb;
    
  public:
    //@{
    //! NTC structures defintion
    typedef std::tuple<double, double, double> vcTup;
    typedef std::map<sKeyPair, vcTup>          vcMap;
    //@}
    
    
    typedef std::map<sKeyPair, QuantileArray> qpMap;
    
  protected:
    
    //@{
    //! Quantile init data
    
    //@{
    //! The desired quantiles.  You can add as many of these 
    //! as you want.
    const double u1 = 0.50, u2 = 0.90, u3 = 0.95, u4 = 0.99;
    const std::vector<double> qs = {u1, u2, u3, u4};
    //@}
    
    //! Main NTC data structures
    qpMap db;
    
    //@{
    //! For debugging
    static unsigned maxrpt;
    static unsigned skprpt;
    key_type caller;
    unsigned report;
    bool inTest;
    
    void VelCrsTest();
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
    static double VelCrsMin;
    
    //! Default <cross section> x <relative velocity> value
    static double VelCrsDef;
    
    //! Return quantile value for a single index
    vcTup VelCrsAvg(sKeyPair indx, double quantile=0.5);
    
    //! Return quantile values for all indices
    vcMap VelCrsAvg(double quantile);
    
    //! Add a single value to the list
    void  VelCrsAdd(sKeyPair indx, const vcTup& val);
    
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
	
	p.second.send();
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
	
	sKeyPair key(p1, p2);
	
	db[key].recv(id);
      }
    }

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

    //! In serialization
    bool inSer;
    
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
      inSer = true;
      
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
      if (!inSer) {
	std::string filename = ".ntcdb." + runtag;
	if ( boost::filesystem::exists(filename) ) {
	  restore_myself(filename);
	}
      }
      inSer = false;
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
  
