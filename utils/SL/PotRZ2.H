#ifndef _PotRZ_H
#define _PotRZ_H

#include "QDHT.H"
#include <cmath>

//! Compute the inverse Hankel transform to get the potential based on
//! Gauss' Law applied to the cylindrical Laplacian.
//!
//! The user supplies a functor containing the surface density.
//! E.g. B&T equations 2.184-2.185.
//!
class PotRZ
{
protected:
  //! Sigma(R) functor instance stored by the constructor
  std::function<double(double)> mass;
  
  //! Compute the Hankel transform of the target density
  class Sk
  {
  protected:

    int N;
    QDHT v;

  public:
    
    //! Null constructor
    Sk() {}

    //! Constructor with R=cutoff, N=number of knots, M=harmonic order,
    //! and a=parameter vector
    Sk(double R, int N, int M) : N(N)
    {
      v = QDHT(M, N, R);
    }
    
    //! Perform the transform
    std::tuple<Eigen::VectorXd, Eigen::VectorXd>
    operator()(std::function<double(double)> dens)
    {
      auto r = v.getR();
      Eigen::VectorXd f(N);
      for (int i=0; i<N; i++) f[i] = dens(r[i]);
      return {v.getK(), v(f)};
    }
  };


private:

  QDHT v;
  Sk   u;

public:

  //! The gravitational potential by Hankel transform
  //!
  //! Constructor 
  PotRZ(double R, int N, int M)
  {
    v = QDHT(M, N, R);
    u = Sk(R, N, M);
  }

  //! Perform the evaluation
  std::tuple<Eigen::VectorXd, Eigen::VectorXd> operator()
  (double z, std::function<double(double)> dens)

  {
    Eigen::VectorXd k, S;
    std::tie(k, S) = u(dens);
    for (int i=0; i<k.size(); i++)
      S[i] = exp(-k[i]*z)/k[i] * S[i] * (-2.0*M_PI);
    return {v.getR(), v(S, false)};
  }

  //! Evalute the forward transform: Sk
  std::tuple<Eigen::VectorXd, Eigen::VectorXd>
  getKT(std::function<double(double)> func)
  {
    return u(func);
  }
};

#endif
