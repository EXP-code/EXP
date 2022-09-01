#ifndef _KMEANS_H
#define _KMEANS_H

#include <algorithm>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <memory>

/*
  This is a simple implementation of Lloyd's k-means algorithm with a
  general distance metric that may be specified using a functor
  interface.
 */
namespace MSSA
{
  namespace KMeans
  {
    /*
      A general concept of a point in an ndim dimensional vector space
      with a cluster ID
    */
    class Point
    {
    public:
      int cid;		   // no default cluster
      int ndim;		   // dimension
      double minDist;	   // default infinite dist to nearest cluster
      std::vector<double> x; // coordinates
      
      //! Null contructor
      Point()
      {
	ndim    = 0;
	cid     = -1;
	minDist = std::numeric_limits<double>::max();
      }
      
      //! Constructor for a particular dimension
      Point(int N)
      {
	ndim    = N;
	cid     = -1;
	minDist = std::numeric_limits<double>::max();
	x.resize(ndim);
      }
      
      //! Construct from data
      Point(const std::vector<double>& p) : x(p)
      {
	ndim    = x.size();
	cid     = -1;
	minDist = std::numeric_limits<double>::max();
      }
      
      //! Access to vector
      double& operator[](int i) { return x[i]; }
      
    };
    
    using Ptr = std::shared_ptr<Point>;
    
    /** Distance functor class
     */
    class KMeansDistance
    {
    public:
      
      virtual double operator()(const std::vector<double>& x,
				const std::vector<double>& y) = 0;
    };
    
    
    class EuclideanDistance : public KMeansDistance
    {
    public:
      double operator()(const std::vector<double>& x,
			const std::vector<double>& y);
    };
    
    class WcorrDistance : public KMeansDistance
    {
    private:
      int numT, numW, Lstar, Kstar;
      
    public:
      /** Constructor
	  @param numT is the size of the time series
	  @param numW is the window length
      */
      WcorrDistance(int numT, int numW) : numT(numT), numW(numW)
      {
	Lstar  = std::min<int>(numT - numW, numW);
	Kstar  = std::max<int>(numT - numW, numW);
      }
      
      //! Compute distance
      double operator()(const std::vector<double>& x,
			const std::vector<double>& y);
    };
    
    class WcorrDistMulti : public KMeansDistance
    {
    private:
      int numT, numW, nchn, Lstar, Kstar;
      
    public:
      /** Constructor
	  @param numT is the size of the time series
	  @param numW is the window length
	  @param nchn is the number of appended channels
      */
      WcorrDistMulti(int numT, int numW, int nchn) :
	numT(numT), numW(numW), nchn(nchn)
      {
	Lstar  = std::min<int>(numT - numW, numW);
	Kstar  = std::max<int>(numT - numW, numW);
      }
      
      //! Compute distance
      double operator()(const std::vector<double>& x,
			const std::vector<double>& y);
    };
    
    
    /**
       Simple implementation of Lloyd's algorithm for arbitrary number of
       dimensions
    */
    class kMeansClustering
    {
    private:
      
      //! The vector of centroids
      std::vector<std::vector<double>> cen;
      
      //! Cache the data dimension for convenience
      int ndim;
      
      //! Copy of the points
      std::vector<Ptr> classes;
      
      //! The last Euclidean difference between centroids
      double total;
      
    public:
      
      //! Constructor
      kMeansClustering(std::vector<Ptr>& points)
      {
	ndim = points[0]->ndim;
	for (auto p : points) classes.push_back(p);
      }
      
      
      /** Perform niter iterations on k means
	  
	  @param dist is the metric distance for grouping
	  @param k is the number of clusters to seed
	  @param s is the stride for center seeding (s>0); random by default
	  @param verbose true prints diagnostic info (false by default)
      */
      void iterate(KMeansDistance& dist, int niter, int k, int s=0,
		   bool verbose=false);
      
      //! Get the centers
      std::vector< std::vector<double> > get_cen() { return cen; }
      
      //! Get the vectors and center ids
      std::vector< std::tuple<std::vector<double>, int> > get_results()
      {
	std::vector< std::tuple<std::vector<double>, int> > ret;
	for (auto p : classes) {
	  ret.push_back({p->x, p->cid});
	}
	return ret;
      }
      
      //! Get the convergence measure: the summed Euclidean distance
      //! between current centroids and previous centroids
      double getTol() { return total; }
      
    };
    
  }
}
  
#endif
  
