#ifndef _DISK_MODELS_H
#define _DISK_MODELS_H

#include <EmpCylSL.h>

//! The usual exponential disk
class Exponential : public EmpCylSL::AxiDisk
{
private:
  double a, h;

public:

  Exponential(double a, double h, double M=1) : a(a), h(h), EmpCylSL::AxiDisk(M) {}

  double operator()(double R, double z)
  {
    double sh = 1.0/cosh(z/h);
    return 0.25*M/(M_PI*a*a*h) * exp(-R/a) * sh * sh;
  }

};

//! Miyamoto-Nagai disk
class MNdisk : public EmpCylSL::AxiDisk
{
private:
  double a, h;

public:

  MNdisk(double a, double h, double M=1) : a(a), h(h), EmpCylSL::AxiDisk(M) {}

  double operator()(double R, double z)
  {
    double Z2 = z*z + h*h;
    double Z  = sqrt(Z2);
    double Q  = a + Z;
    return 0.25*h*h*M/M_PI*(a*R*R + (a + 3.0*Z)*Q*Q)/( pow(R*R + Q*Q, 2.5) * Z * Z2 );
  }

};

#endif
