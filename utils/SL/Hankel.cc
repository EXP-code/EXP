#define _USE_MATH_DEFINES // using M_PI for pi

#include <iostream>
#include <cmath>

#include "Hankel.H"

// Constructor
HankelTransform::HankelTransform(double h, double nu, int N)
{
  // Use mapped formula by default
  mapped = true;

  // Check the step
  if ( h >= 0.0) {
    this->h = h;
  } else {
    std::ostringstream sout;
    sout << "h (" << h << ") must be positive";
    throw std::runtime_error(sout.str());
  }

  // Check the Bessel order
  if ( nu >= 0.0) {
    this->nu = nu;
  } else {
    std::ostringstream sout;
    sout << "nu (" << nu << ") must be positive";
    throw std::runtime_error(sout.str());
  }

  // Check the number of knots
  if ( N >= 1) {
    this->N = N;
  } else {
    std::ostringstream sout;
    sout << "N (" << N << ") must be greater than zero";
    throw std::runtime_error(sout.str());
  }

  // Sets maximum number of nodes to about 2^15:
  //
  const int maxN = 32769;
  
  // Imports zeros of the Bessel function. Initializing this way speeds up calls
  //
  try {
    boost::math::cyl_bessel_j_zero(this->nu, 1, maxN, std::back_inserter(jn_zeros0));
    for (size_t i = 0; i < maxN; i++) {
      zeros.push_back( jn_zeros0[i] );
      xi.push_back( zeros[i]/M_PI );

      // The functions gsl_sf_bessel_Jn and gsl_sf_bessel_Yn return
      // the result of the Bessel functions of the first and second
      // kinds respectively
      //
      Jp1.push_back( boost::math::cyl_bessel_j(nu+1.,M_PI*xi[i]) );
      w.push_back( boost::math::cyl_neumann(nu,M_PI*xi[i])/Jp1[i] );  
    }
  }
  catch (std::exception& ex) {
    std::cout << "Thrown exception " << ex.what() << std::endl;
  }

};

// From Equation 5.1 in Ogata's paper
double HankelTransform::get_psi(double t)
{
  return t*tanh(0.5*M_PI * sinh(t));
};

// From Equation 5.2 in Ogata's paper
double HankelTransform::get_psip(double t)
{
  return M_PI*t*(-pow(tanh( 0.5*M_PI*sinh(t)), 2.0) + 1.0)*0.5*cosh(t) + tanh(0.5*M_PI*sinh(t));
};


// Transformed Ogata quadrature sum. Equation 5.2 in Ogata's paper
//
double HankelTransform::ogata_transformed
(std::function<double (double) > f, double q, double h)
{
  double nu = this->nu;
  int N = this->N;

  double val = 0;

  try {
    for (size_t i = 0; i < (unsigned)N; i++) {
      double knots = M_PI/h*get_psi( h*xi[i] );
      double Jnu = boost::math::cyl_bessel_j(nu,knots);
      double temp =  get_psip( h*xi[i] ), psip;

      // Sanity check
      if (std::isnan(temp)) {
        psip = 1.0;
      }
      else
      {
        psip = temp;
      };

      double F = fk_trans(knots, f, q);
      val += M_PI * w[i]* F * Jnu * psip;
    }
  }
  catch (std::exception& ex)
  {
    std::cout << "Thrown exception " << ex.what() << std::endl;
  }

  return val;
};

// Untransformed, linear Ogata quadrature sum. Equation 1.1 in Ogata's
// paper.
//
double HankelTransform::ogata_linear
(std::function<double (double) > f, double q, double h){
  double nu = this->nu;
  int N = this->N;

  double val = 0;

  try
  {
    for (size_t i = 0; i < (unsigned)N; i++) {
      double knots = xi[i]*h;
      double F = fk_trans(knots, f, q)*boost::math::cyl_bessel_j(nu, knots);
      val += h * w[i] * F;
    }
  }
  catch(std::exception& ex)
  {
    std::cout << "Thrown exception " << ex.what() << std::endl;
  }

  return val;
};

// Compute the Ogata quadrature
double HankelTransform::operator()(std::function<double (double) > g,  double q)
{
  double result = 0.0;

  // The default: Ogata formula transformed to the infinite interval
  if (mapped){
    result = ogata_transformed(g, q, h);
  }
  else {
    result = ogata_linear(g, q, h);
  }

  return result;
};
