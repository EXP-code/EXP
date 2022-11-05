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
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <tuple>
#include <cmath>

// MDW classes
//
#include <numerical.H>
#include <gaussQ.H>
#include <EmpCylSL.H>
#include <DiskModels.H>

#include <config_exp.h>

double  plgndr(int l, int m, double x);
double dplgndr(int l, int m, double x);

/** Compute the axisymmetric force of a three-dimensional disk density
    using the multipole expansion

    @param model is the target model for potential computation
    @param rmin is the minimum radius for potential evaluation
    @param rmax is the maximum radius for potential evaluation
    @param ascl is the scale radius for maping r to the unit interval.
    This is only used if ascl>0, otherwise it is linear if rmin=0 or
    logarithmic if rmin>0.
    @param lmax is the maximum order of the spherical harmonics
    @param numr is the number of grid points in the radial table
    @param nint is the number of quadrature knots for the x=cos(theta) integral
    @parm progress set to true prints a progress bar

    \note{Mathematics}

    Using the standard multipole formula
    \f[
    \Phi(\mathbf{r}) = -4\pi G \sum_{lm} Y_{lm}(\mathbf{\Omega})
    \left[
    \frac{1}{r^{l+1}}\int_0^r da a^{l+2}\rho_{lm}(a) +
    r^l\int_r^{\infty} da \frac{\rho_{lm}(a)}{a^{l-1}}
    \right]
    \f]
    where
    \f[
    \rho_{lm} = \int d\mathbf{\Omega}Y_{lm}(\mathbf{\Omega}) \rho(\mathbf{r}).
    \f]

    The forces follow by analytic differentiation. Specifically, the
    forces follows by differentiation of \f$\Phi(\mathbf{r})\f$
    with respect to \f$r\f$ and \f$\theta\f$ where:
    \f[
    F_R = -\frac{\partial\Phi}{\partial r} \frac{R}{r} - 
    \frac{\partial\Phi}{\partial\theta} \frac{z}{r^2}
    \f]
    and
    \f[
    F_z = -\frac{\partial\Phi}{\partial r} \frac{z}{r} +
    \frac{\partial\Phi}{\partial\theta} \frac{R}{r^2}
    \f]
*/    
class DiskEval
{
private:
  EmpCylSL::AxiDiskPtr model;
  double xmin, xmax, dx, ascl;
  int    lmax, numr;
  bool   logr, xscl;

  static const std::string cachefile;

  //! Rho(l, r), Term1(l, r), Term2(l, r)
  //std::vector< std::vector<double> > rho, T1, T2;

  //! Rho(l, m, r), Term1(l, m, r), Term2(l, m, r)
  std::vector< std::vector< std::vector<double> > > rho, T1, T2;

  //! Ylm angular factors
  // CAUTION: this will fail at l>150 owing to double precision issues in exp call.
  double Nlm(int l, int m)
  {
    return sqrt( (2.0*l+1.0)/(4.0*M_PI) ) * 
      exp(0.5*(lgamma(1.0+l-m) - lgamma(1.0+l+m)));
  }

  //! Ylm evaluation
  // CAUTION: this will fail at l>150 owing to double precision issues.
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

  double x_to_r(double x)
  {
    return ascl * x/(1.0 - x);
  }

  double r_to_x(double r)
  {
    return r/(r + ascl);
  }

  double dr_to_dx(double x)
  {
    return ascl/(1.0 - x)/(1.0 - x);
  }

  bool read_cache();
  void write_cache();

  const unsigned hmagic = 0x57ee1ed;
  
public:

  //! Constructor
  DiskEval(EmpCylSL::AxiDiskPtr model, double rmin, double rmax,
	   double ascl, int lmax, int numr, int nint, bool progress=true, int mmax=1, int nump=1, bool cache=true);

  //! Evaluation where the return tuple is potential, dPhi/dR, dPhi/dz, dPhi/dphi
  std::tuple<double, double, double, double> operator()(double R, double z, double phi=0);

};

#endif

