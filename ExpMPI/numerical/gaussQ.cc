
#pragma implementation

#include <math.h>
#include <string>
#include <Vector.h>
#include <gaussQ.h>

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

  w.setsize(1, n);
  r.setsize(1, n);

  Hermite( n, alpha, r.array(0,n-1), w.array(0,n-1) );

  w *= exp(lgamma(0.5+0.5*alpha));
}

LaguQuad::LaguQuad(int N, double ALPHA)
{
  FunctionID = "LaguQuad";

  n = N;
  alpha = ALPHA;

  w.setsize(1, n);
  r.setsize(1, n);

  Laguerre( n, alpha, r.array(0,n-1), w.array(0,n-1) );

  w *= exp(lgamma(1.0+alpha));
}

JacoQuad::JacoQuad(int N, double ALPHA, double BETA)
{
  FunctionID = "JacoQuad";

  n = N;
  alpha = ALPHA;
  beta = BETA;

  w.setsize(1, n);
  r.setsize(1, n);

  Jacobi( n, alpha, beta, r.array(0,n-1), w.array(0,n-1) );

  w *= exp(lgamma(1.0+alpha) + lgamma(1.0+beta) - lgamma(2.0+alpha+beta));
}

