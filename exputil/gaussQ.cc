
#include <cstdlib>
#include <cmath>
#include <string>

#include "gaussQ.H"

extern "C" {
  void Hermite(int n, double alpha, double abscis[], double weight[]);
  void Jacobi(int n, double alpha, double beta, 
	      double abscis[], double weight[]);
  void Laguerre(int n, double alpha, double abscis[], double weight[]);
};

HermQuad::HermQuad(int N, double ALPHA)
{
  FunctionID = "HermQuad";

  n = N;
  alpha = ALPHA;

  w.resize(n);
  r.resize(n);

  Hermite( n, alpha, r.data(), w.data() );

  for (auto & v : w) v *= exp(lgamma(0.5+0.5*alpha));
}

LaguQuad::LaguQuad(int N, double ALPHA)
{
  FunctionID = "LaguQuad";

  n = N;
  alpha = ALPHA;

  w.resize(n);
  r.resize(n);

  Laguerre( n, alpha, r.data(), w.data() );

  for (auto & v : w) v *= exp(lgamma(1.0+alpha));
}

JacoQuad::JacoQuad(int N, double ALPHA, double BETA)
{
  FunctionID = "JacoQuad";

  n = N;
  alpha = ALPHA;
  beta = BETA;

  w.resize(n);
  r.resize(n);

  Jacobi( n, alpha, beta, r.data(), w.data() );

  for (auto & v : w )
    v *= exp(lgamma(1.0+alpha) + lgamma(1.0+beta) - lgamma(2.0+alpha+beta));
}

