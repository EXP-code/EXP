#ifndef _Quantile_H
#define _Quantile_H

#include <cmath>
#include <vector>
#include <algorithm>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

namespace NTC {

  class Quantile
  {
    friend class boost::serialization::access;

  private:
    const int ssize = 5;

    std::vector<double> dn, npos, hgt;
    std::vector<int> pos;

    double p;
    unsigned long num;

    bool full;

    //! Piecewise parabolic computation
    double P2(double d,
	      double qp1, double q, double qm1, 
	      double np1, double n, double nm1);

    //! Update the heights and positions according to the P2 algorithm
    void update();
  
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(dn);
      ar & BOOST_SERIALIZATION_NVP(npos);
      ar & BOOST_SERIALIZATION_NVP(hgt);
      ar & BOOST_SERIALIZATION_NVP(pos);
      ar & BOOST_SERIALIZATION_NVP(p);
      ar & BOOST_SERIALIZATION_NVP(num);
      ar & BOOST_SERIALIZATION_NVP(full);
    }

  public:

    //! Reset quantile value
    void reset(double P);

    //! Null constructor
    Quantile() { reset(0.5); }

    //! Constructs a single quantile object
    Quantile(double p) { reset(p); }

    //! Add another datum
    void operator()(double item);
    
    //! Return current quantile
    double operator()();

    //! Return data count so far
    unsigned long count() { return num; }

    //@{
    //! MPI syncronize
    
    //! Node sends data to root
    void send();

    //! Root intializes itself from node's data
    void recv(int id);
    //@}
    
  };

}

#endif
