#include <iostream>
#include <cmath>

#include "QDHT.H"

bool QDHT::debug = false;

// Constructor
QDHT::QDHT(int nu, int N, double R) : nu(nu), N(N), R(R)
{
  // Check the Bessel order
  if ( nu >= 0.0) {
    this->nu = nu;
  } else {
    std::ostringstream sout;
    sout << "nu (" << nu << ") must be positive";
    throw std::runtime_error(sout.str());
  }

  // Check the number of knots
  if ( N >= 1) {
    this->N = N;
  } else {
    std::ostringstream sout;
    sout << "N (" << N << ") must be greater than zero";
    throw std::runtime_error(sout.str());
  }

  // Sets maximum number of nodes to about 2^15:
  //
  const int maxN = 32769;
  
  // Imports zeros of the Bessel function. Initializing this way speeds up calls
  //
  try {
    boost::math::cyl_bessel_j_zero(this->nu, 1, maxN, std::back_inserter(zeros));
  }
  catch (std::exception& ex) {
    std::cout << "Thrown exception " << ex.what() << std::endl;
  }

  // Assign the total bandwidth: R*V
  S = zeros[N];

  // Deduce V from S and R
  V = S/R;

  // Assign the workspace
  r.resize(N);
  k.resize(N);
  T.resize(N, N);
  Jp.resize(N);

  // Get the Bessel function of order nu+1 evaluated at the roots
  for (int i=0; i<N; i++) {
    Jp[i] = boost::math::cyl_bessel_j(nu+1.0, zeros[i]);
  }

  // Precompute the radial and spatial frequency vectors and the
  // transform matrix
  for (int i=0; i<N; i++) {
    r(i) = zeros[i]/V;
    k(i) = zeros[i]/R;

    for (int j=i; j<N; j++) {

      T(i, j) = 2.0/S *
	boost::math::cyl_bessel_j(nu, zeros[i]*zeros[j]/S) / (Jp[i]*Jp[j]);
      if (i != j) T(j, i) = T(i, j);
    }
  }

  if (debug) check();

};

