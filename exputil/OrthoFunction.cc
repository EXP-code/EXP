#include <OrthoFunction.H>

OrthoFunction::OrthoFunction
(int nmax, DensFunc W, double rmin, double rmax, double scl, int dof, bool seg) :
  nmax(nmax), W(W), rmin(rmin), rmax(rmax), scale(scl), dof(dof), segment(seg)
{
  // Quadrature knots/weights
  lq = std::make_shared<LegeQuad>(knots);

  // Set up for inner product computation
  xmin = r_to_x(rmin);
  xmax = r_to_x(rmax);
  dx   = xmax - xmin;
  
  // Compute the recursion coefficients
  generate();
}

void OrthoFunction::generate()
{
  alph.resize(nmax+1);
  beta.resize(nmax+1);
  norm.resize(nmax+1);

  beta[0] = norm[0] = scalar_prod(0, 0);
  alph[0] = scalar_prod(0, 1)/norm[0];

  for (int i=1; i<=nmax; i++) {
    norm[i] = scalar_prod(i, 0);
    alph[i] = scalar_prod(i, 1)/norm[i];
    beta[i] = norm[i]/norm[i-1];
  }
}

double OrthoFunction::scalar_prod(int n, int m)
{
  double ans = 0.0;

  for (int i=0; i<knots; i++) {
    double x = xmin + dx*lq->knot(i);
    double r = x_to_r(x);
    double y = 2.0*r/scale; if (segment) y = x;
    Eigen::VectorXd p = poly_eval(y, n) * W(r);
    ans += dx*lq->weight(i)*d_r_to_x(x)*pow(r, dof-1) *
      p(n) * p(n) * pow(y, m);
  }

  return ans;
}

Eigen::VectorXd OrthoFunction::poly_eval(double t, int n)
{
  Eigen::VectorXd ret(n+1);
  ret[0] = 1.0;
  if (n) {
    ret[1] = (t - alph[0])*ret[0];
    for (int j=1; j<n; j++) {
      ret[j+1] = (t - alph[j])*ret[j] - beta[j]*ret[j-1];
    }
  }

  return ret;
}

Eigen::MatrixXd OrthoFunction::testOrtho()
{
  Eigen::MatrixXd ret(nmax+1, nmax+1);
  ret.setZero();

  for (int i=0; i<knots; i++) {
    // Get abscissa and weight
    double x = xmin + dx*lq->knot(i);
    double r = x_to_r(x);
    double f = dx*lq->weight(i)*d_r_to_x(x)*pow(r, dof-1);
    double y = 2.0*r/scale; if (segment) y = x;

    // Evaluate polynomial
    Eigen::VectorXd p = poly_eval(y, nmax) * W(r);

    // Compute scalar products
    for (int n1=0; n1<=nmax; n1++) {
      for (int n2=0; n2<=nmax; n2++) {
	ret(n1, n2) += f * p(n1) * p(n2) / sqrt(norm[n1]*norm[n2]);
      }
    }
  }

  return ret;
}

Eigen::VectorXd OrthoFunction::operator()(double r)
{
  // Enforce bounds
  r = std::max<double>(rmin, std::min<double>(rmax, r));

  // Weight
  double w = W(r);

  // Choose system;
  double y = 2.0*r/scale;
  if (segment) y = r_to_x(r);

  // Generate normalized orthogonal functions
  auto p = poly_eval(y, nmax);
  for (int i=0; i<=nmax; i++) p(i) *= w/sqrt(norm[i]);

  return p;
}
