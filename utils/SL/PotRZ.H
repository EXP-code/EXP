#ifndef _PotRZ_H
#define _PotRZ_H

#include "Hankel.H"

//! Compute the inverse Hankel transform to get the potential.  The
//! user supplies a functor containing the (mass per cylindrical
//! radius>/(2*pi).  E.g. surface density time radius.
class PotRZ
{
public:
  //! The dM(R)/dR functor
  std::function<double(double)> mass;
  
protected:
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
  //! h1 and h2 are the steps for the inverse and Sk transform
  //! N1 and N2 are the knot sizes for the inverse and Sk
  //! M is the harmonic order
  //! a is the model parameter vector
  PotRZ(double h1, double h2, int N1, int N2, int M,
	std::function<double(double)> func) : mass(func)
  {
    v = HankelTransform(h1, (double)M, N1);
    u = Sk(h2, N2, M, this);
  }

  //! Perform the evaluation
  double operator()(double R, double z)
  {
    return v(std::bind(&PotRZ::func, this, std::placeholders::_1, z), R);
  }
};

#endif
