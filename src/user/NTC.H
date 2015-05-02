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
  typedef std::map<sKeyPair, unsigned>       numMap;
  //@}

protected:

  //@{
  //! NTC structures
  dqMap     VelCrsList;
  vcMap     VelCrsSum;
  numMap    VelCrsNum;
  //@}
  
  //@{
  //! For debugging
  static unsigned maxrpt;
  key_type caller;
  unsigned report;
  void VelCrsTest();
  void setKey(const key_type& k) { caller = k; }
  //@}

  // Maximum number of cached relative velocities
  static size_t VelCrsSZ; 

public:

  //! Null constructor
  NTCitem() { report = 0; }

  //! Copy constructor
  NTCitem(boost::shared_ptr<NTCitem> p) 
  {
    VelCrsList = p->VelCrsList;	
    VelCrsSum  = p->VelCrsSum;	
    VelCrsNum  = p->VelCrsNum;    
    report     = 0;
  }

  //! Threshold <cross section> x <relative velocity> value
  static double VelCrsMin;

  //! Default <cross section> x <relative velocity> value
  static double VelCrsDef;

  //@{
  //! Maximum relative velocity average (<0 if no entries)
  vcTup VelCrsAvg(sKeyPair indx);
  vcMap VelCrsAvg();
  //@}
  
  //@{
  //! Add maximum to list
  void VelCrsAdd(vcMap& vals);
  void VelCrsAdd(sKeyPair indx, const vcTup& val);
  //@}
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
