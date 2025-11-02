#ifndef _KDdensity_H_
#define _KDdensity_H_

#include <functional>
#include <cstdint>

#include <Eigen/Eigen>

#include "ParticleReader.H"
#include "permutation.H"
#include "localmpi.H"
#include "KDtree.H"

namespace Utility
{
  class KDdensity
  {
  protected:
    using point3 = KDtree::point <double, 3>;
    using tree3  = KDtree::kdtree<double, 3>;

    //! Point list
    std::vector<point3> points;

    //! Density array
    std::vector<double> KDdens;

    //! Accumulated mass of phase space
    double KDmass = 0.0;

    //! Number of points per ball
    int Ndens = 32;

    //! The KD tree instance
    std::shared_ptr<tree3> kdtree_;

    //! The point index to particle ID map
    std::map<unsigned long, unsigned> indx;
    
  public:
    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    //! Create KD tree from a particle reader
    KDdensity(PR::PRptr reader, int Ndens=32);

    //! Create KD tree from a arrays of masses and positions
    KDdensity(const Eigen::VectorXd& m, const RowMatrixXd& d, int Ndens=32);

    //! Get density for particle index
    double getDensity(unsigned long index);

  };
}

#endif
