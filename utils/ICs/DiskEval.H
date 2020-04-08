// Do multipole evaluation of 3-dimensional axisymmetric disk density
//
#ifndef _DiskEval_H
#define _DiskEval_H

                                // C++/STL headers
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <tuple>
#include <cmath>

// MDW classes
//
#include <numerical.h>
#include <gaussQ.h>
#include <EmpCylSL.h>
#include <DiskModels.H>

#include "config.h"

double  plgndr(int l, int m, double x);
double dplgndr(int l, int m, double x);

class DiskEval
{
private:
  EmpCylSL::AxiDiskPtr model;
  double rmin, rmax, dr;
  int    lmax, numr;
  bool   logr;

  //! Rho(l, r), Term1(l, r), Term2(l, r)
  std::vector< std::vector<double> > rho, T1, T2;

  //! Ylm angular factors
  double Nlm(int l, int m)
  {
    return sqrt( (2.0*l+1.0)/(4.0*M_PI) ) * 
      exp(0.5*(lgamma(1.0+l-m) - lgamma(1.0+l+m)));
  }

  //! Ylm evaluation
  double Ylm(int l, int m, double cosx)
  {
    int M = abs(m);
    double plm = plgndr(l, M, cosx);
    if (std::isnan(plm) or std::isinf(plm))
      std::cout << "Failure in pldndr at l=" << l
		<< " m=" << m << " cosx=" << cosx << std::endl;
    return Nlm(l, M) * plm * pow(-1.0, M);
  }

  //! Partial derivative of Ylm in theta
  double Zlm(int l, int m, double cosx)
  {
    if (l==0 or fabs(cosx)>=1.0) return 0.0;

    int M = abs(m);
    double dplm = dplgndr(l, M, cosx);
    if (std::isnan(dplm) or std::isinf(dplm))
      std::cout << "Failure in pldndr at l=" << l
		<< " m=" << m << " cosx=" << cosx << std::endl;

    return -Nlm(l, M) * dplm * pow(-1.0, M) * sqrt(1.0 - cosx*cosx);
  }

public:

  //! Constructor
  DiskEval(EmpCylSL::AxiDiskPtr model, double rmin, double rmax,
	   int lmax, int numr, int nint, bool progress=true);

  //! Evaluation where the return tuple is potential, dphi/dR, dphi/dz
  std::tuple<double, double, double> operator()(double R, double z);

};

#endif
