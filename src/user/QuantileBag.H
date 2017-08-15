#ifndef _QuantileBag_H
#define _QuantileBag_H

#include "Quantile.H"

namespace NTC {

  /**
     @brief Container of Quantile objects

     Maintains an array 5-point quantiles and an equally spaced histogram
  */
  class QuantileBag
  {
    friend class boost::serialization::access;

  protected:
    //@{
    //! Quantile init data

    std::map<double, Quantile> quant;
    Quantile hist;
    unsigned int M;
    
    //@{
    //! The desired quantiles.  You can add as many of these 
    //! as you want.
    const int Nequal = 10;
    const double u1 = 0.01, u2 = 0.05, u3 = 0.95, u4 = 0.99;
    const std::vector<double> qs = {u1, u2, u3, u4};
    //@}
    
    //! Count instances for debugging only
    static unsigned instance;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(quant);
      ar & BOOST_SERIALIZATION_NVP(hist);
      ar & BOOST_SERIALIZATION_NVP(M);
    }

  public:

    //! Null constructor
    QuantileBag();

    //! Only for debugging instances
    ~QuantileBag();

    //! Copy constructor
    QuantileBag(const QuantileBag&);

    //! Copy operator
    QuantileBag &operator=(const QuantileBag&);

    //! Add another datum
    //! @param item is the data to be added
    void add(double item);
    
    //! Get the quantile with value p
    double operator()(double p);
    
    //! Return data count so far
    unsigned long count() { return M; }

    //! Get the value probability value for \f$x\f$, i.e. \f$P^{-1}(x)\f$
    double inverse(double x);

    //! Get minimum value so far
    double xmin();

    //! Get maximum value so far
    double xmax();

    //@{
    //! MPI syncronize
    
    //! Node sends data to root
    void send();

    //! Root intializes itself from node's data
    void recv(int id);
    //@}
    
    //! For debugging: dump accumulated data points
    void dump(std::ostream& out);

    //! Return if quantile estimator is fully populated
    bool full() { return(M >= 5); }

    //! Number of live instances (for debugging)
    static unsigned live() { return instance; }
  };

}

#endif
