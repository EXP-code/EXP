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
  bool   logr, xscl, prenorm=true;

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


  //! Standard recursion
  //! CAUTION: this will fail at l>150 owing to double precision issues.
  double Ylm_standard(int l, int m, double cosx)
  {
    int M = abs(m);
    double plm = plgndr(l, M, cosx);
    if (std::isnan(plm) or std::isinf(plm)) // Check for failure in plgndr
      std::cout << "Failure in plgndr at l=" << l
		<< " m=" << m << " cosx=" << cosx << std::endl;
    return Nlm(l, M) * plm * pow(-1.0, M);
  }

  //! Partial derivative of Ylm in theta
  //! Standard version
  double Zlm_standard(int l, int m, double cosx)
  {
    if (l==0 or fabs(cosx)>=1.0) return 0.0;

    int M = abs(m);
    double dplm = dplgndr(l, M, cosx);
    if (std::isnan(dplm) or std::isinf(dplm))
      std::cout << "Failure in dplgndr at l=" << l
		<< " m=" << m << " cosx=" << cosx << std::endl; // Correct function name is dplgndr

    return -Nlm(l, M) * dplm * pow(-1.0, M) * sqrt(1.0 - cosx*cosx);
  }


  //! Prenormalized Ylm evaluation
  double Ylm_prenorm(int L, int m, double x)
  {
    int M = abs(m);
  
    // Initial value precompile
    constexpr double ylm00 = 0.25/M_2_SQRTPI;
  
    // Assign initial value
    double pmm = ylm00;
  
    // Recurrence to get l=m=M
    if (M>0) {
      double somx2 = std::sqrt( (1.0 - x)*(1.0 + x) );
      for (int l=1; l<=M; l++) {
	pmm = -std::sqrt(1.0 + 0.5/l)*somx2*pmm;
      }
    }
    
    // Even or odd M?
    double facM = 1.0;
    if (M & 0x1) facM = -1.0;
  
    // We are done if L equals M
    if (L == M) {
      return pmm * facM;
    }
    else {
      // L=M+1
      double pmmp1 = std::sqrt(1.0 + 2.0*(M+1))*x*pmm, pll = 0.0;
      if (L == M+1) {
	return pmmp1 * facM;
      }
      // L>M+1 recursion
      else {
	for (int l=M+2; l<=L; l++) {
	  // Recursion constants
	  double pfac  = (2.0*l + 1.0)/(2.0*l - 3.0)/(l*l - M*M);
	  double alpha = std::sqrt( pfac * (4*(l-1)*(l-1) - 1 ) );
	  double beta  = std::sqrt( pfac * ((l-1)*(l-1) - M*M ) );
	  // Compute next term
	  pll   = alpha*x*pmmp1 - beta*pmm;
	  pmm   = pmmp1;
	  pmmp1 = pll;
	}
	return pll * facM;
      }
    }
  }
  
  //! Derivative of Ylm prenormalized
  double Zlm_prenorm(int L, int m, double x)
  {
    int M = abs(m);
    double Ylmp = Ylm_prenorm(L+1, M, x), Ylm = Ylm_prenorm(L, M, x);
    double fac = -1.0/std::sqrt(1.0 - x*x);
    double fac2 = 2.0*L + 1.0, fac3 = 2.0*L + 3.0;
    return fac*( Ylm*x*(1.0 + L) - Ylmp*std::sqrt( fac2/fac3*(fac2 + L*L - M*M) ) );
  }

  //! Ylm evaluation
  double Ylm(int l, int m, double cosx)
  {
    if (prenorm)
      return Ylm_prenorm(l, m, cosx);
    else
      return Ylm_standard(l, m, cosx);
  }

  //! Ylm derivative
  double Zlm(int l, int m, double cosx)
  {
    if (prenorm)
      return Zlm_prenorm(l, m, cosx);
    else
      return Zlm_standard(l, m, cosx);
  }

  //! Coordinate scaling: scaled to physical
  double x_to_r(double x)
  {
    return ascl * x/(1.0 - x);
  }

  //! Coordinate scaling: physical to scaled
  double r_to_x(double r)
  {
    return r/(r + ascl);
  }

  // Coordinate scaling: Jacobian
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

  //! Set prenorm (true) standard (false) evaluation
  void setPrenorm(bool b) { prenorm = b; }

};

#endif

