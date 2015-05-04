#ifndef _NTC_H
#define _NTC_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <deque>
#include <map>

#include <boost/shared_ptr.hpp>

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
  friend class NTCdb;

public:
  //@{
  //! NTC structures defintion
  typedef std::tuple<double, double, double> vcTup;
  typedef std::deque<vcTup>                  dqTup;
  typedef std::map<sKeyPair, dqTup>          dqMap;
  typedef std::map<sKeyPair, vcTup>          vcMap;
  //@}

  //! Compute a running quantile using a relaxation scheme
  class Quantile
  {
  private:

    double eta;
    vcTup  q;
    double p;
    unsigned long c;

  public:

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

    //! Add value
    void add(const vcTup& v)
    {
      std::get<0>(q) += eta * (std::copysign(1.0, std::get<0>(v) - std::get<0>(q)) + 2.0*p - 1.0);
      std::get<1>(q) += eta * (std::copysign(1.0, std::get<1>(v) - std::get<1>(q)) + 2.0*p - 1.0);
      std::get<2>(q) += eta * (std::copysign(1.0, std::get<2>(v) - std::get<2>(q)) + 2.0*p - 1.0);
      c++;
    }

    //! Get current value
    vcTup get()
    {
      return q;
    }

    //! Get number of samples
    unsigned long count()
    {
      return c;
    }
  };

  typedef boost::shared_ptr<Quantile>  qPtr;

  class QuantileArray
  {
  private:
    std::vector<double>    pp;
    std::vector<vcTup>     qq;
    std::map<double, qPtr> tt;

  public:
    //! Null constructor
    QuantileArray() {}

    //! Construct an array of quantile pointers
    QuantileArray
    (const std::vector<std::tuple<double, vcTup, double> >& e)
    {
      initialize(e);
    }

    //! Initialize the array
    void initialize
    (const std::vector<std::tuple<double, vcTup, double> >& e)
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
  };
  
  typedef std::map<sKeyPair, QuantileArray> qpMap;

protected:

  //@{
  //! Quantile init list
  typedef std::tuple<double, vcTup, double> d3Tup;
  const std::set<double> qs = {0.05, 0.50, 0.95};
  const double eta = 0.1, vi = 1.0;
  const std::vector<d3Tup>
  ic = { d3Tup(eta, vcTup(vi, vi, vi), 0.05),
	 d3Tup(eta, vcTup(vi, vi, vi), 0.50), 
	 d3Tup(eta, vcTup(vi, vi, vi), 0.95) };
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

  // Maximum number of cached relative velocities
  static size_t VelCrsSZ; 

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
  vcMap VelCrsAvg(double quantile=0.5);
  
  //! Add multiple data values to the list
  void VelCrsAdd(dqMap& vals);

  //! Add a single value to the list
  void VelCrsAdd(sKeyPair indx, const vcTup& val);
};

typedef boost::shared_ptr<NTCitem> NTCptr;
typedef std::map<key_type, NTCptr> NTCdata;

//! Provide a nice interface for the NTCdata std::map
class NTCdb
{
protected:

  //! The std::map
  NTCdata data;

public:

  //! Return shared pointer to NTCitem instance for cell key k
  NTCptr operator[](const key_type& k);

};

#endif
