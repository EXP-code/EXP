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

  };

  typedef boost::shared_ptr<Quantile> qPtr;

  class QuantileArray
  {
  public:

    typedef std::tuple<double, vcTup, double> d3;
    typedef std::vector<d3> vd3;

  private:
    std::map<double, qPtr> tt;

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
  };
  
  typedef std::map<sKeyPair, QuantileArray> qpMap;

protected:

  //@{
  //! Quantile init data

  //@{
  //! The 3 desired quantiles
  const double u1 = 0.25, u2 = 0.50, u3 = 0.95;
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
  
  //! Add multiple data values to the list
  void VelCrsAdd(const dqMap& vals);

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
