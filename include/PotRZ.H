#ifndef _PotRZ_H
#define _PotRZ_H

#include <QDHT.H>
#include <cmath>

//! Compute the inverse Hankel transform using the quasi-discrete
//! Hankel transform to get the potential based on Gauss' Law applied
//! to the cylindrical Laplacian.  The user supplies a functor
//! containing the surface density.  E.g. B&T equations 2.184-2.185.
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

  QDHT v, vm, vp;
  Sk   u;
  int  M;

public:

  //! The gravitational potential by Hankel transform
  //!
  //! Constructor 
  PotRZ(double R, int N, int M) : M(M)
  {
    v = QDHT(M, N, R);
    u = Sk(R, N, M);
    if (M>0) vm = QDHT(M-1, N, R);
    vp = QDHT(M+1, N, R);
  }

  //! Which field
  enum class Field { potential, rforce, zforce };

  //! Perform the evaluation
  std::tuple<Eigen::VectorXd, Eigen::VectorXd> operator()
  (double z, std::function<double(double)> dens, Field f=Field::potential)

  {
    Eigen::VectorXd k, S, Sm, Sp;
    std::tie(k, S) = u(dens);

    if (f == Field::rforce) { Sm.resize(k.size()); Sp.resize(k.size()); }

    // Prepare the transform for the desired output field
    //
    for (int i=0; i<k.size(); i++) {
      switch (f) {
      case Field::potential:
	S[i] = -exp(-k[i]*z)/k[i] * S[i] * 2.0*M_PI;
	break;
      case Field::zforce:
	S[i] = -exp(-k[i]*z) * S[i] * 2.0*M_PI;
	break;
      case Field::rforce:
	if (M>0) Sm[i] = -exp(-k[i]*z) * S[i] * 2.0*M_PI;
	Sp[i] = exp(-k[i]*z) * S[i] * 2.0*M_PI;
	break;
      }
    }

    // Compute the inverse transform
    //
    if (f==Field::potential or f==Field::zforce) {
      // Potential or zforce
      return {v.getR(), v(S, false)};
    } else {
      // Derivative recursion
      if (M==0) return {v.getR(), -vp(Sp, false)};
      return {v.getR(), (vm(Sm, false) - vp(Sp, false))*0.5};
    }
  }

  //! Perform the evaluation at a particular r
  double operator()
  (double r, double z, std::function<double(double)> dens, Field f=Field::potential)
  {
    Eigen::VectorXd k, S, Sm, Sp;
    std::tie(k, S) = u(dens);

    if (f == Field::rforce) { Sm.resize(k.size()); Sp.resize(k.size()); }

    // Prepare the transform for the desired output field
    //
    for (int i=0; i<k.size(); i++) {
      switch (f) {
      case Field::potential:
	S[i] = -exp(-k[i]*z)/k[i] * S[i] * 2.0*M_PI;
	break;
      case Field::zforce:
	S[i] = -exp(-k[i]*z) * S[i] * 2.0*M_PI;
	break;
      case Field::rforce:
	if (M>0) Sm[i] = -exp(-k[i]*z) * S[i] * 2.0*M_PI;
	Sp[i] = exp(-k[i]*z) * S[i] * 2.0*M_PI;
	break;
      }
    }

    // Perform the inverse
    //
    if (f==Field::potential or f==Field::zforce) {
      // Potential or zforce
      return v(r, S);
    } else {
      // Derivative recursion
      if (M==0) return -vp(r, Sp);
      return (vm(r, Sm) - vp(r, Sp))*0.5;
    }
  }

  //! Evalute the forward transform: Sk
  std::tuple<Eigen::VectorXd, Eigen::VectorXd>
  getKT(std::function<double(double)> func)
  {
    return u(func);
  }
};

#endif
