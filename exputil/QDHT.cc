#include <iostream>
#include <cmath>

#include "QDHT.H"

bool QDHT::debug = false;

Eigen::VectorXd bessjz(int n, int m);

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

  // Imports zeros of the Bessel function. Initializing this way speeds up calls
  //
  try {
    zeros  = bessjz(nu, N+1);
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
    Jp[i] = EXPmath::cyl_bessel_j(nu+1.0, zeros[i]);
  }

  // Precompute the radial and spatial frequency vectors and the
  // transform matrix
  for (int i=0; i<N; i++) {
    r(i) = zeros[i]/V;
    k(i) = zeros[i]/R;

    for (int j=i; j<N; j++) {

      T(i, j) = 2.0/S *
	EXPmath::cyl_bessel_j(nu, zeros[i]*zeros[j]/S) / (Jp[i]*Jp[j]);
      if (i != j) T(j, i) = T(i, j);
    }
  }

  if (debug) check();

}

Eigen::VectorXd QDHT::operator()(Eigen::VectorXd& v, bool forward)
{
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


double QDHT::operator()(double r, Eigen::VectorXd& v)
{
  double ret = 0.0;
  for (int i=0; i<N; i++) {
    ret +=  2.0/(R*R*Jp[i]*Jp[i])*v[i]*EXPmath::cyl_bessel_j(nu, zeros[i]*r/R);
  }
  return ret;
}

void QDHT::check()
{
  double det = T.determinant();
  std::cout << "QDHT: solution quality=" << std::fabs(det-1.0)
	    << " det=" << det << std::endl;
  
  // Perform the SVD
  //
  auto begin = std::chrono::steady_clock::now();
  Eigen::BDCSVD<Eigen::MatrixXd>
    svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
  auto S = svd.singularValues();
  auto end = std::chrono::steady_clock::now();
  auto tim = std::chrono::duration<double>(end-begin).count();
  
  int good=0, bad=0;
  double worst = 0.0, det1 = 1.0;
  for (size_t i=0; i<S.size(); i++) {
    det1 *= S(i);
    double test = fabs(S(i) - 1.0);
    if (test < 1.0e-8) {
      good++;
    } else {
      if (test>worst) worst = test;
      bad++;
    }
  }
  
  std::cout << "QDHT: singular values: good=" << good << " bad=" << bad
	    << " worst=" << worst << " det=" << det1
	    << " time(SVD)=" << tim << " sec"
	    << std::endl << std::endl;
}
