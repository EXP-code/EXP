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

#include <boost/shared_ptr.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>

#include <boost/filesystem.hpp>

#include "pHOT_types.H"
#include "Species.H"

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
// Subtract two tuples (non const version)
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

  //! Compute a running quantile using a relaxation scheme
  class Quantile
  {
    friend class boost::serialization::access;
    
  private:

    double eta;
    vcTup  q;
    double p;
    unsigned long c;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(eta);
      ar & BOOST_SERIALIZATION_NVP(q);
      ar & BOOST_SERIALIZATION_NVP(p);
      ar & BOOST_SERIALIZATION_NVP(c);
    }

  public:

    //! Null constructor
    Quantile() : eta(0.0), p(0.0), c(0) {}

    //! Constructor
    Quantile(double eta, const vcTup& q, double p) : 
      eta(eta), q(q), p(p), c(0) {}

    //! Copy constructor
    Quantile(const Quantile *s)
    {
      eta = s->eta;
      q   = s->q;
      p   = s->p;
      c   = s->c;
    }

    /** Add value: the main algorithm
	
	The expectation of [sgn(v-q)] is zero at the median (p=0.5).

	The expectation of [sgn(v-q) + 2p - 1] is zero at the
	p-quantile.  That is, if q is the p-quantile for the
	distribution of v, then v-q is smaller than 0 with probability
	p.  Or:

	<[]> = <sgn(v-q)> + 2p - 1 = (-1)*p + 1*(1 - p) + 2p - 1 = 0.

	Thus, q is larger than the quantile p, (v-q) will be negative
	with probability > p and [] will be more often negative than
	positive (and vice versa). In this way, each value is
	"relaxed" toward the expectation by the fraction eta of the
	current median estimate.
     */
    void add(const vcTup& v)
    {
      // Do each tuple entry separately
      //
      std::get<0>(q) += eta * fabs(std::get<0>(q)) *
	(std::copysign(1.0, std::get<0>(v) - std::get<0>(q)) + 2.0*p - 1.0);

      std::get<1>(q) += eta * fabs(std::get<1>(q)) *
	(std::copysign(1.0, std::get<1>(v) - std::get<1>(q)) + 2.0*p - 1.0);

      std::get<2>(q) += eta * fabs(std::get<2>(q)) *
	(std::copysign(1.0, std::get<2>(v) - std::get<2>(q)) + 2.0*p - 1.0);

      // For diagnostic output
      c++;
    }

    //! Get current value
    const vcTup& get() { return q; }

    //! Get number of samples
    const unsigned long count() { return c; }

    //@{
    //! MPI syncronize

    //! Node sends data to root
    void send()
    {
      MPI_Send(&eta,              1, MPI_DOUBLE, 0, 1101, MPI_COMM_WORLD);
      MPI_Send(&(std::get<0>(q)), 1, MPI_DOUBLE, 0, 1102, MPI_COMM_WORLD);
      MPI_Send(&(std::get<1>(q)), 1, MPI_DOUBLE, 0, 1103, MPI_COMM_WORLD);
      MPI_Send(&(std::get<2>(q)), 1, MPI_DOUBLE, 0, 1104, MPI_COMM_WORLD);
      MPI_Send(&p,                1, MPI_DOUBLE, 0, 1105, MPI_COMM_WORLD);
      MPI_Send(&c,                1, MPI_UNSIGNED_LONG, 0, 1106, MPI_COMM_WORLD);
    }

    // Root intializes itself from node's data
    void recv(int id)
    {
      MPI_Recv(&eta,              1, MPI_DOUBLE, id, 1101, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&(std::get<0>(q)), 1, MPI_DOUBLE, id, 1102, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&(std::get<1>(q)), 1, MPI_DOUBLE, id, 1103, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&(std::get<2>(q)), 1, MPI_DOUBLE, id, 1104, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&p,                1, MPI_DOUBLE, id, 1105, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&c, 1, MPI_UNSIGNED_LONG, id, 1106, MPI_COMM_WORLD, 
	       MPI_STATUS_IGNORE);
    }
    //@}

  };

  typedef boost::shared_ptr<Quantile> qPtr;

  class QuantileArray
  {
    friend class boost::serialization::access;

  public:

    typedef std::tuple<double, vcTup, double> d3;
    typedef std::vector<d3> vd3;

  private:
    std::map<double, qPtr> tt;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(tt);
      if (Archive::is_loading::value) {
	std::cout << ".................." << std::endl;
	std::cout << " Loading [" << myid << "]" << std::endl;
	std::cout << ".................." << std::endl;
	std::cout << "Size=" << tt.size() << std::endl;
	for (auto v : tt) {
	  std::cout << "v=" << v.first    << std::endl;
	}
	std::cout << ".................." << std::endl;
      }

      if (Archive::is_saving::value) {
	std::cout << ".................." << std::endl;
	std::cout << " Saving [" << myid << "]" << std::endl;
	std::cout << ".................." << std::endl;
	std::cout << "Size=" << tt.size() << std::endl;
	for (auto v : tt) {
	  std::cout << "v=" << v.first    << std::endl;
	}
	std::cout << ".................." << std::endl;
      }
    }

  public:
    //! Null constructor
    QuantileArray() {}

    //! Construct an array of quantile pointers
    QuantileArray(const vd3& e)
    {
      initialize(e);
    }

    //! Initialize the array
    void initialize(const vd3& e)
    {
      tt.erase(tt.begin(), tt.end());

      for (auto v : e) {
	double eta = std::get<0>(v);
	vcTup  q   = std::get<1>(v);
	double p   = std::get<2>(v);
	tt[p] = qPtr(new Quantile(eta, q, p));
      }
    }
    
    void add(const vcTup& v) 
    {
      for (auto &q : tt) q.second->add(v);
    }

    qPtr operator[](double p)
    {
      return tt[p];
    }

    size_t size()
    {
      return tt.size();
    }

    //! Node sends to root
    void send()
    {
      // Sends the array of Quantiles
      //
      unsigned siz = tt.size();
      MPI_Send(&siz, 1, MPI_UNSIGNED, 0, 1130, MPI_COMM_WORLD);
      for (auto &p : tt) {
	double d = p.first;
	MPI_Send(&d, 1, MPI_DOUBLE, 0, 1131, MPI_COMM_WORLD);
	p.second->send();
      }
    }
      
    //! Root initializes itself from node's data
    void recv(int id) 
    {
      unsigned siz;
      MPI_Recv(&siz, 1, MPI_UNSIGNED, id, 1130, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      // Receives the array of quantiles from a node
      //
      for (unsigned j=0; j<siz; j++) {
	double d;
	MPI_Recv(&d, 1, MPI_DOUBLE, id, 1131, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);

	qPtr q(new Quantile);
	q->recv(id);
	tt[d] = q;
      }

      if (tt.size()==0) {
	std::cout << "Root received zero tt size from Node " << id 
		  << std::endl;
      }
    }

  }; // end: QuantileArray
  
  typedef std::map<sKeyPair, QuantileArray> qpMap;

protected:

  //@{
  //! Quantile init data

  //@{
  //! The 3 desired quantiles
  const double u1 = 0.50, u2 = 0.90, u3 = 0.95;
  const std::set<double> qs = {u1, u2, u3};
  //@}

  //@{
  //! The precision and initial tuple values
  const double et = 0.1, vi = 1.0;
  const vcTup v3 = vcTup(vi, vi, vi);
  //@}

  //@{
  //! The constant intialization vector
  typedef std::tuple<double, vcTup, double> d3;
  typedef std::vector<d3> vd3;
  const vd3 ic = {d3(et, v3, u1), d3(et, v3, u2), d3(et, v3, u3)};
  //@}
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

public:

  //! Constructor
  NTCitem() 
  { 
    report = 0; 
    inTest = false; 
  }

  //! Copy constructor
  NTCitem(boost::shared_ptr<NTCitem> p) 
  {
    db     = p->db;
    report = 0;
    inTest = false;
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
  void VelCrsAdd(sKeyPair indx, const vcTup& val);

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
      
      if (p.second.size()==0) {
	std::cout << "Sending zero db size in Node " << myid 
		  << " for key=" << p.first << std::endl;
      }
    }
    
  }

  // Root node initializes from node
  //
  void recv(int id)
  {
    MPI_Recv(&caller, 1, MPI_UNSIGNED_LONG, id, 1110, MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
      
    unsigned siz;

    MPI_Recv(&siz, 1, MPI_UNSIGNED,      id, 1111, MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);

    for (unsigned k=0; k<siz; k++) {

      speciesKey p1, p2;

      MPI_Recv(&p1.first,  1, MPI_UNSIGNED_SHORT, id, 1112, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&p1.second, 1, MPI_UNSIGNED_SHORT, id, 1113, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&p2.first,  1, MPI_UNSIGNED_SHORT, id, 1114, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      MPI_Recv(&p2.second, 1, MPI_UNSIGNED_SHORT, id, 1115, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      sKeyPair key(p1, p2);
      
      db[key].recv(id);
      
      if (db[key].size()==0) {
	std::cout << "Received zero db size in Node " << myid 
		  << " for key=" << key << std::endl;
      }
    }
  }

};

typedef boost::shared_ptr<NTCitem> NTCptr;
typedef std::map<key_type, NTCptr> NTCdata;

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


//! Provide a nice interface for the NTCdata std::map
class NTCdb
{
  friend class boost::serialization::access;

protected:

  //! The std::map
  NTCdata data;

  //! Last serialization step
  int last_step;

  //! Next serialization step
  int next_step;
  
  //! Save interval
  const int intvl = 1;

  //! In serialization
  bool inSer;

  //! Persistent copy
  NTCdata pdata;

  //! Get copy from all nodes to save
  void sync() 
  {
    if (myid==0) pdata = data;

    unsigned siz;
    unsigned long k;

    for (int i=1; i<numprocs; i++) {
      if (i==myid) {
	siz = data.size();
	MPI_Send(&siz, 1, MPI_UNSIGNED, 0, 1110, MPI_COMM_WORLD);
	for (auto &p : data) {
	  k = p.first;
	  MPI_Send(&k, 1, MPI_UNSIGNED_LONG, 0, 1111, MPI_COMM_WORLD);
	  p.second->send();
	}

      }

      if (0==myid) {
	MPI_Recv(&siz, 1, MPI_UNSIGNED, i, 1110, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	for (unsigned j=0; j<siz; j++) {
	  MPI_Recv(&k, 1, MPI_UNSIGNED_LONG, i, 1111, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

	  // Get the instance from the other node
	  //
	  NTCptr q(new NTCitem);
	  q->recv(i);	  

	  key_type key(k);
	  pdata[key] = q;
	}
      }

    }
  }

  //! Save contents to a file
  void save_myself(const std::string& filename)
  {
    // Synchronize
    sync();

    // make an archive
    std::ofstream ofs(filename.c_str());
    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(pdata);
    oa << BOOST_SERIALIZATION_NVP(last_step);
    oa << BOOST_SERIALIZATION_NVP(next_step);
  }

  //! Restore contents from a file
  void restore_myself(const std::string& filename)
  {
    inSer = true;

    // open the archive
    std::ifstream ifs(filename.c_str());
    boost::archive::xml_iarchive ia(ifs);

    // restore the data from the archive
    ia >> BOOST_SERIALIZATION_NVP(data);
    ia >> BOOST_SERIALIZATION_NVP(last_step);
    ia >> BOOST_SERIALIZATION_NVP(next_step);
  }

public:

  //! Constructor
  NTCdb()
  {
    last_step = this_step;
    next_step = last_step + intvl;
    if (!inSer) {
      std::string filename = ".ntcdb." + runtag;
      if ( boost::filesystem::exists(filename) ) {
	restore_myself(filename);
      }
    }
    inSer = false;
  }

  //! Return shared pointer to NTCitem instance for cell key k
  NTCptr operator[](const key_type& k);

  //! Update the serialization file
  void update()
  {
    if (this_step >= next_step) {
      std::string filename = ".ntcdb." + runtag;
      save_myself(filename);
      next_step = this_step + intvl;
    }
  }
  
};

#endif
