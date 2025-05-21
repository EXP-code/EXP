#ifndef _OrthoFunction_H_
#define _OrthoFunction_H_

#include <functional>
#include <memory>

#include <Eigen/Eigen>

#include <gaussQ.H>

//! Compute orthogonal functions whose lowest-order member is the
//! background density using the Stieltjes procedure for generating
//! 3-term recursion
class OrthoFunction
{
public:
  
  using DensFunc = std::function<double(double)>;

private:

  //! Sqrt root weight function
  DensFunc W;

  //! Polynomial order
  int nmax;

  //! Mapping scale
  double scale;

  //! Degrees of freedom
  int dof;

  //! Radial boundaries
  double rmin, rmax;

  //! Mapped radial boundaries and range
  double xmin, xmax, dx;

  //! Recursion coefficients and norm
  Eigen::VectorXd alph, beta, norm;

  //! Legendre knots
  int knots = 400;

  //! Gaussian quadrature
  std::shared_ptr<LegeQuad> lq;

  //@{
  //! Coordinate scaling
  
  double r_to_x(double r)
  {
    double ra = r/scale;
    return (ra - 1.0)/(ra + 1.0);
  }

  double x_to_r(double x)
  {
    return scale * (1.0 + x)/(1.0 - x);
  }

  double d_r_to_x(double x)
  {
    double mx = 1.0 - x;
    return 2.0 * scale / ( mx * mx );
  }
  //@}

  //! Inner product polynomial of order n with moment m
  double scalar_prod(int n, int m);

  //! Polynomial evaluation at t for orders j=0,...,n
  Eigen::VectorXd poly_eval(double t, int n);

  //! Generate the recursion for the orthogonal polynomials
  void generate();

public:

  /** Constructor
  @param norder  Order of the orthogonal polynomial
  @param W       Density weight function
  @param rmin    Inner radial bound
  @param rmax    Outer radial bound
  @param scale   Mapping scale
  @param dof     Degrees of freedom (2 or 3)
  */
  OrthoFunction(int norder, DensFunc W, double rmin, double rmax, double scale, int dof);

  //! Destructor
  virtual ~OrthoFunction() {}

  //! Evaluate orthogonal functions at r
  //! @param r  Radial coordinate
  //! @return   A vector of the orthogonal functions
  Eigen::VectorXd operator()(double r);

  //! Reset Legendre knots and weights for inner product
  //! @param N  Number of Legendre knots
  void setKnots(int N)
  {
    knots = N;
    lq = std::make_shared<LegeQuad>(N);
    generate();
  }

  //! Test orthogonality
  Eigen::MatrixXd testOrtho();

  //! Dump the orthogonal function table
  //! @param filename  Filename for output
  void dumpOrtho(const std::string& filename);
};


#endif
