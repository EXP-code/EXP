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

     This version implements both the single quantile algorithm from
     Box 1 and the multiple quantile algorithm from Box 2.

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
    static const size_t ssize = 5;
    
    //! Set true for debugging and sanity checking MPI
    static const bool DBG_VERBOSE = false;

    //! Count instances for debugging only
    static unsigned instance;

    //! Saved data points
    std::vector<double> q;
    
    //! Marker values
    std::vector<double> dn;

    //! Position values for markers
    std::vector<double> np;

    //! Current indices of saved data points
    std::vector<int>    n;
    
    //! Number of datums so far
    int M;

    //! Number of markers
    int N;

    //@{
    //! For storage management and updating
    void    init();
    double *extend(int count);
    void    update();
    double  parabolic( int i, int d );
    double  linear   ( int i, int d );
    //@}

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(q);
      ar & BOOST_SERIALIZATION_NVP(dn);
      ar & BOOST_SERIALIZATION_NVP(np);
      ar & BOOST_SERIALIZATION_NVP(n);
      ar & BOOST_SERIALIZATION_NVP(M);
      ar & BOOST_SERIALIZATION_NVP(N);

      if (DBG_VERBOSE) {
	std::string stat("loading");
	if (Archive::is_saving::value) stat = "saving";
	
	if (q.size() != ssize) {
	  int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	  std::cout << "[" << myid
		    << "] crazy q size=" << q.size()
		    << " on " << stat << std::endl;
	  dn.resize(ssize);
	}

	if (dn.size() != ssize) {
	  int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	  std::cout << "[" << myid
		    << "] crazy dn size=" << dn.size()
		    << " on " << stat << std::endl;
	  dn.resize(ssize);
	}

	if (np.size() != ssize) {
	  int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	  std::cout << "[" << myid
		    << "[ crazy np size=" << np.size()
		    << " on " << stat << std::endl;
	  np.resize(ssize);
	}

	if (n.size() != ssize) {
	  int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	  std::cout << "[" << myid
		    << "] crazy n size=" << n.size()
		    << " on " << stat << std::endl;
	  n.resize(ssize);
	}
      }
    }

  public:

    //! Add a new quantile value
    //! @param P the quantile probability
    void newQ(double P);

    //! Make a histogram with \f$n\f$ bins
    void histogram(int n);

    //! Null constructor
    Quantile();

    //! Constructs a single quantile object
    //! @param P the quantile probability
    Quantile(double p);

    //! Copy constructor
    //! @param p another Quantile instance
    Quantile(const Quantile& p);

    //! Currently, only for debugging
    ~Quantile();

    //! Add another datum
    //! @param item is the data to be added
    void add(double item);
    
    //! Get the quantile
    double operator()();
    
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

    //@{
    //! Get internal data count and target (for debugging only)
    int datums() { return M; }
    int target() { return N; }
    //@}

    //! Return if quantile estimator is fully populated
    bool full() { return M >= N; }


    //! Number of live instances (for debugging)
    static unsigned live() { return instance; }
  };

}

#endif
