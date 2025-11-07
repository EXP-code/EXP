#ifndef _QDHT_H
#define _QDHT_H

#include <functional>
#include <algorithm>
#include <chrono>
#include <vector>
#include <cmath>

#include <Eigen/Eigen>

#include "EXPmath.H"		// For C++17 special function workarounds

//! Hankel transform based on the QDHT algorithm
//!
//! See:
//! "Computation of quasi-discrete Hankel transforms of the integer
//! order for propagating optical wave fields" Manuel Guizar-Sicairos
//! and Julio C. Guitierrez-Vega J. Opt. Soc. Am. A 21 (1) 53-58
//! (2004)

class QDHT
{
private:
  //! nu is Bessel function order
  double nu;

  //! N is number of knots
  int N;

  //! Will contain the zeros of the Bessel function of order nu
  Eigen::VectorXd zeros;

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
  //! For checking unitarity (default: false)
  static bool debug;

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
  Eigen::VectorXd operator()(Eigen::VectorXd& v, bool forward=true);

  //! Inverse Hankel transform at r
  double operator()(double r, Eigen::VectorXd& v);

  //! Compute the input coordinates (e.g. radii)
  const Eigen::VectorXd& getR() { return r; }

  //! Compute the transform coordinates (e.g. spatial frequencies)
  const Eigen::VectorXd& getK() { return k; }

  //! Check solution
  void check();
};

#endif
// END: _QDHT_H
