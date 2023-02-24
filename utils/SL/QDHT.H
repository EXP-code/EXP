#ifndef _QDHT_H
#define _QDHT_H

#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>

#include <Eigen/Eigen>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/tools/minima.hpp>

//! Hankel transform based on the QDHT algorithm
//!
class QDHT
{
private:
  //! nu is Bessel function order
  double nu;

  //! N is number of knots
  int N;

  //! Will contain the zeros of the Bessel function of order nu
  std::vector<double> zeros;

  //! Dimension 1 scale
  double R;

  //! Total bandwidth
  double S;

  //! Dimension 2 scale
  double V;

  //! R and K vectors
  Eigen::VectorXd r, k;

  //! Jp vector
  Eigen::VectorXd Jp;

  //! Transform matrix
  Eigen::MatrixXd T;

public:
  //! Null constructor
  QDHT() {}

  //! Constructor
  //! @param M is the harmonic order
  //! @param N is the number of knots
  //! @param R is the radius at which the input function is zero
  QDHT(int M, int N, double R);

  //! Destructor
  ~QDHT() {}

  //! Compute the forward or backward Hankel transform
  Eigen::VectorXd operator()(Eigen::VectorXd& v, bool forward=true) {
    if (forward) {
      Eigen::VectorXd F = v.array() / Jp.array() * R;
      auto G = T * F;
      return G.array() * Jp.array() / V;
    } else {
      Eigen::VectorXd G = v.array() / Jp.array() * V;
      auto F = T * G;
      return F.array() * Jp.array() / R;
    }
  }

  //! Compute the input coordinates (e.g. radii)
  const Eigen::VectorXd& getR() { return r; }

  //! Compute the transform coordinates (e.g. spatial frequencies)
  const Eigen::VectorXd& getK() { return k; }
};

#endif
// END: _QDHT_H
