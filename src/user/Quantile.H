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

#include <mpi.h>

namespace NTC {

  /**
     @brief Estimate a quantile using the Jain & Chlamtac P^2 algorithm

     The \f$P^2\f$ algorithm estimates a quantile without storing
     samples. Five points, or "markers" as the authors called them are
     maintained and sequentially updated based on new data.  The
     values or "heights" of these markers are the minimum and the
     maximum of the samples and the current estimates of the
     \f$(p/2)\f$-, \f$p\f$- and \f$(1+p)/2\f$-quantiles. Their
     positions are equal to the number of samples that are smaller or
     equal to the markers. Each time a new samples is recorded, the
     positions of the markers are updated and if necessary their
     heights are adjusted using a piecewise- parabolic formula.

     For further details, see:

     R. Jain and I. Chlamtac, The P^2 algorithm for dynamic
     calculation of quantiles and histograms without storing
     observations, Communications of the ACM, Volume 28 (October),
     Number 10, 1985, p. 1076-1085.
    */
  class Quantile
  {
    friend class boost::serialization::access;

  private:
    const size_t ssize = 5;

    static unsigned instance;

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
    //! @param P the quantile probability
    void reset(double P);

    //! Null constructor.  Assumes quantile probabilty of 0.5
    //! (i.e. median) by default.
    Quantile() { 
      reset(0.5); 
      instance++; 
    }

    //! Constructs a single quantile object
    //! @param P the quantile probability
    Quantile(double p) { 
      reset(p); 
      instance++;
    }

    //! Copy constructor
    //! @param p another Quantile instance
    Quantile(const Quantile& p);

    //! Currently, only for debugging
    ~Quantile();

    //! Add another datum
    //! @param item is the data to be added
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
