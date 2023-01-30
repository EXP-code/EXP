#ifndef _CauchyPV_H
#define _CauchyPV_H

// 
// <summary> Header for Cauchy Principal Value integration using 
//           Legendre polynomials </summary>
//
//  Could be generalized (easily) for other sets of polies.  E.g. will do
//  Jacobi polies as soon as references are obtained from Interlibrary Loan
//  Cauchy principal value integrals computed using quadrature rules
//  following Paget and Elliott (Numer. Math. 19, 373).
//
//  Integral performed is:
//
//  <srcblock>
//
//         1
//         /            1
//  I(L) = | dx f(x) -------
//         /          L - x
//        -1
//
//  </srcblock>
//
//  To evaluate integral:
//
//  <srcblock>
//
//         b
//         /            1
//  I(M) = | dy g(y) -------
//         /          M - y
//         a
//
//
//  </srcblock>
//
//  use the transformation:
//
//  <srcblock>
//
//
//         2         b + a
//   L = -----  M -  -----
//       b - a       b - a
//
//  </srcblock>
//
//  with
//
//  <srcblock>
//
//             b - a       b + a
//   f(x) = g( ----- x  +  ----- )
//               2           2
//
//  </srcblock>
//
//
//

#include <vector>
#include <Eigen/Eigen>

class PVQuad : public Legendre
{
private:
  int n;
  Eigen::VectorXd eign;
  std::vector<Eigen::VectorXd> poly;
  Eigen::VectorXd chrs;
  Eigen::VectorXd nrml;
  Eigen::VectorXd coefs;

  void get_eigenvalues(void);
  void get_orthopolies(void);
  Eigen::VectorXd& get_coefs(Eigen::VectorXd&);
  double get_integral(double);

  double q0(double x) {return log( (1.0+x)/(1.0-x) );}

public:

  // Create quadrature with N knots
  PVQuad(int N=10);

  // After constructtion, use root() to access knots and
  double root(int i) {return eign[i];}

  // use mu() get access weights
  double mu(int i) { return chrs[i];}
				// Scalar product norm
  double norm(int i) { return nrml[i];}
				// Orthogonal polynomial
  Eigen::VectorXd& ovec(int i) {return poly[i];}

				// Coefficients for input function
  Eigen::VectorXd& return_coefs(std::function<double(double)> func);
  Eigen::VectorXd& return_coefs(Eigen::VectorXd& input) {return get_coefs(input);}

				// Cauchy Principal Value integral
				// for given pole

  double return_pv(double lambda, std::function<double(double)> func) {
    return_coefs(func); return get_integral(lambda);}
    
  double return_pv(double lambda, Eigen::VectorXd& input) {
    return_coefs(input); return get_integral(lambda);}

  double return_pv(double lambda) { return get_integral(lambda);}

};

#endif
