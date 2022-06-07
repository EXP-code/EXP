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


#define BINARY_ARCHIVE		// Comment for xml archive

#ifdef BINARY_ARCHIVE
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#else
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#endif

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>

#include <boost/version.hpp>

#include "pHOT_types.H"
#include "QuantileBag.H"
#include "Species.H"
#include "global.H"		// For simulation time: tnow




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


namespace NTC
{
  //! Interaction type with interaction elements
  typedef std::tuple<unsigned short, speciesKey, speciesKey> T;

  //! Default null
  extern speciesKey const nullKey;

  //! Default electron
  extern speciesKey const electron;
  
  //! Default proton
  extern speciesKey const proton;

  //! Default interaction key
  extern T const single;

  //! Container for interaction list to add UI to data types
  template <class U>
  class Interact
  {

  public:

    //! Array-like data implemented by a map
    std::map<T, U> v;

    //! Single datum
    U d;

    //! Constructor (for initialization)
    Interact() : d(0) {}

    //! Is map empty?
    bool operator!() { return v.size() == 0; }

    //! Get datum from the map
    U& operator[](T t) { return v[t]; }

    //! Get single datum
    U& operator()() { return d; }

  };

  class NTCitem
  {
    friend class boost::serialization::access;
    friend class NTCdb;
    
  public:

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

    class uqMap : public std::map<T, QuantileBag>
    {
      friend class boost::serialization::access;

    public:

      uqMap() {}

      uqMap(const QuantileBag& q) 
      {
	(*this)[single] = QuantileBag();
	(*this)[single] = q;
      }

      uqMap(T t, const QuantileBag& q)
      {
	(*this)[t] = QuantileBag();
	(*this)[t] = q;
      }

      //! Invoked for debugging data values before/after serialization
      void dump(std::ostream& out)
      {
	for (auto v : *this) {
	  out << std::string(72, '-') << std::endl
	      << std::string(20, '=') << " Interact " << v.first << std::endl;
	  v.second.dump(out);
	}
      }
      
    private:

      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
      {
	typedef std::map<T, QuantileBag> Tup3QuantileBagMap;
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Tup3QuantileBagMap);
      }
    };

    typedef std::map<sKeyPair, udMap> vcMap;
    typedef std::map<sKeyPair, uqMap> qpMap;
    //@}
    
  protected:
    
    //! Main NTC data structures
    qpMap db;
    
    //! Default time stamp
    static double def_time;

    //! Time stamp
    double ts;

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
      ar & BOOST_SERIALIZATION_NVP(ts);
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

    class Error : public std::exception 
    {
    public:
      Error() {}
      virtual ~Error() throw() {}
      virtual const char *what() const throw()
      {
	return "NTCitem: no CrsVel available";
      }
    };
    
  public:
    
    //! Constructor
    NTCitem() 
    { 
      ts     = def_time;
      report = 0; 
      inTest = false;

      instance++;
    }
    
    //! Copy constructor
    NTCitem(const NTCitem& p) 
    {
      ts     = p.ts;
      db.insert(p.db.begin(), p.db.end());

      report = 0;
      inTest = false;

      instance++;
    }
    
    //! Assignment operator
    NTCitem &operator=(const NTCitem& p)
    {
      ts     = p.ts;
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

    //! Invoked for debugging data values before/after serialization
    void dump(std::ostream& out);

    //! Node sends to root
    void send()
    {
      static unsigned long dnn = 0;

      dnn++;
      MPI_Send(&dnn, 1, MPI_UNSIGNED_LONG, 0, 1108, MPI_COMM_WORLD);


      unsigned long clr = caller;

      unsigned siz      = db.size();
      
      MPI_Send(&ts,  1, MPI_DOUBLE,        0, 1109, MPI_COMM_WORLD);
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
	speciesKey tk;

	for (auto &q : p.second) {
	  MPI_Send(&(u=std::get<0>(q.first)),         1, MPI_UNSIGNED_SHORT, 0, 1117, MPI_COMM_WORLD);

	  MPI_Send(&(u=std::get<1>(q.first).first),   1, MPI_UNSIGNED_SHORT, 0, 1118, MPI_COMM_WORLD);
	  tk = std::get<1>(q.first);
	  MPI_Send(&tk.first,                         1, MPI_UNSIGNED_SHORT, 0, 1119, MPI_COMM_WORLD);
	  MPI_Send(&tk.second,                        1, MPI_UNSIGNED_SHORT, 0, 1120, MPI_COMM_WORLD);

	  MPI_Send(&(u=std::get<2>(q.first).first),   1, MPI_UNSIGNED_SHORT, 0, 1121, MPI_COMM_WORLD);
	  tk = std::get<2>(q.first);
	  MPI_Send(&tk.first,                         1, MPI_UNSIGNED_SHORT, 0, 1122, MPI_COMM_WORLD);
	  MPI_Send(&tk.second,                        1, MPI_UNSIGNED_SHORT, 0, 1123, MPI_COMM_WORLD);

	  q.second.send();
	}
      }
      
      // Compare resent tag for sync check
      unsigned long dnn1;
      MPI_Recv(&dnn1, 1, MPI_UNSIGNED_LONG, 0, 1124, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (dnn != dnn1) {
	std::cout << "Node " << myid << " synchronization error, got="
		  << dnn1 << " but expected " << dnn << std::endl;
      }

    }
    
    // Root node initializes from node
    //
    void recv(int id)
    {
      unsigned long dnn = 0;
      MPI_Recv(&dnn,    1, MPI_UNSIGNED_LONG, id, 1108, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&ts,     1, MPI_DOUBLE,        id, 1109, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
	  speciesKey tk;
	  T t;
	  MPI_Recv(&(get<0>(t)),        1, MPI_UNSIGNED_SHORT, id, 1117, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  MPI_Recv(&(get<1>(t).first),  1, MPI_UNSIGNED_SHORT, id, 1118, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&tk.first,           1, MPI_UNSIGNED_SHORT, id, 1119, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&tk.second,          1, MPI_UNSIGNED_SHORT, id, 1120, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  std::get<1>(t) = tk;

	  MPI_Recv(&(get<2>(t).first),  1, MPI_UNSIGNED_SHORT, id, 1121, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&tk.first,           1, MPI_UNSIGNED_SHORT, id, 1122, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&tk.second,          1, MPI_UNSIGNED_SHORT, id, 1123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  std::get<2>(t) = tk;

	  db[key][t].recv(id);

	} // db2 size loop

      } // db1 size loop

      // Resend tag for sync check
      MPI_Send(&dnn, 1, MPI_UNSIGNED_LONG, id, 1124, MPI_COMM_WORLD);

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
    
    //! Copy data to ascii file
    static bool debug_output;

    //! Ascii version number
    static unsigned debug_count;

    //! Get copies from all nodes to save
    void sync() 
    {
      // Insert root's nodes into pdata
      if (myid==0) pdata.insert(data.begin(), data.end());
      
      unsigned siz;

      unsigned long k;
      
      for (int i=1; i<numprocs; i++) {
	if (i==myid) {
	  siz = data.size();
	  MPI_Send(&siz, 1, MPI_UNSIGNED, 0, 1110, MPI_COMM_WORLD);
	  for (auto &p : data) {
	    k = p.first;
	    MPI_Send(&k,  1, MPI_UNSIGNED_LONG, 0, 1111, MPI_COMM_WORLD);
	    p.second.send();
	  }
	  
	}
	
	if (0==myid) {
	  MPI_Recv(&siz, 1, MPI_UNSIGNED, i, 1110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for (unsigned j=0; j<siz; j++) {
	    MPI_Recv(&k,  1, MPI_UNSIGNED_LONG, i, 1111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    
	    // Get the instance from the other node
	    //

	    NTCitem tmp;
	    tmp.recv(i);

	    NTCdata::iterator it = pdata.find(k);
	    if (it == pdata.end()) {
	      pdata[k] = tmp;
	    } else {
	      if (tmp.ts > it->second.ts) it->second = tmp;
	    }
	  }
	}
	// Need this barrier to prevent the MPI buffer cache from
	// getting insanely large
	//
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    
    //! Save contents to a file
    void save_myself(const std::string& filename);
    
    //! Restore contents from a file
    void restore_myself(const std::string& filename);
    
    //! Invoked for debugging data values before/after serialization
    void dump();

  public:
    
    //! Chatty output for debugging
    static bool chatty;

    //! Save interval
    static unsigned intvl;
    
    //! Maximum age for uncalled cell entries in number of time steps
    static unsigned maxAge;

    //! Constructor
    NTCdb();
    
    //! Return NTCitem instance for cell key k
    NTCitem& operator[](const key_type& k);
    
    //! Update the serialization file
    void update()
    {
      if (tnow >= next_time) {
	// Set <current time>
	curr_time = tnow;
	// Set time for next save
	next_time = curr_time + dtime*(static_cast<double>(intvl)-eps);
	std::string filename = ".ntcdb." + runtag;
	save_myself(filename);
      }
    }
    
    //! Serialize to finish
    void finish()
    {
      // Set current time but not next save time
      curr_time = tnow;
      std::string filename = ".ntcdb." + runtag;
      save_myself(filename);
    }

    //! Clean up unused cell entries
    void ageOut()
    {
      double delT = dtime * maxAge;
      for (auto it=data.begin(); it!=data.end();) {
	if (tnow - it->second.ts > delT) it = data.erase(it);
	else ++it;
      }
    }

  };

  //! Maintains a maximum value as a mapped object
  class maxDouble
  {
  private:
    double value;

  public:
    //! Constructor and assign initial value
    maxDouble(double init=0.0) : value(init) {}

    //! Maximum operator
    void operator[](double val)
    {
      if (val>value) value = val;
    }

    //! Get current value
    double& operator()()
    {
      return value;
    }

  };

  typedef Interact<maxDouble> InteractD;

} // End: namespace NTC

#endif
  
