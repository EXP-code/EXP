#ifndef _PotRZ_H
#define _PotRZ_H

#include <Hankel.H>

//! Compute the inverse Hankel transform to get the potential based on
//! Gauss' Law applied to the cylindrical Laplacian.
//!
//! The user supplies a functor containing the <mass per cylindrical
//! radius>/(2*pi) which is the surface density times the cylindrical
//! radius.  E.g. B&T equations 2.184-2.185.
//!
class PotRZ
{
protected:
  //! The R*Sigma(R) functor instance stored by the constructor
  std::function<double(double)> mass;
  
  //! Compute the Hankel transform of the target density
  class Sk
  {
  protected:

    HankelTransform v;
    PotRZ *ptr;

  public:
    
    //! Null constructor
    Sk() {}

    //! Constructor with h=step, N=number of knots, M=harmonic order,
    //! and a=parameter vector
    Sk(double h, int N, int M, PotRZ *ptr) : ptr(ptr)
    {
      v = HankelTransform(h, (double)M, N);
    }
    
    //! Perform the transform
    double operator()(double k)
    {
      return -2.0*M_PI*v(ptr->mass, k);
    }
  };


private:

  HankelTransform v;
  Sk u;

  double func(double k, double z)
  {
    return u(k) * exp(-k*z);
  }

public:

  //! The gravitational potential by Hankel transform
  //!
  //! @param h1 is the step count for the inverse transform
  //! @param h2 is the steps cout for the Sk transform
  //! @param N1 is the knot number for the inverse
  //! @param N2 is the knot number for Sk
  //! @param M is the harmonic order
  //! @param a is the model parameter vector
  PotRZ(double h1, double h2, int N1, int N2, int M,
	std::function<double(double)> func) : mass(func)
  {
    v = HankelTransform(h1, (double)M, N1);
    u = Sk(h2, N2, M, this);
  }

  //! Change the functor without recomputing the knots
  void setFunction(std::function<double(double)> func) { mass = func; }

  //! Perform the evaluation
  double operator()(double R, double z)
  {
    return v(std::bind(&PotRZ::func, this, std::placeholders::_1, z), R);
  }
};

#endif
