#include <limits>
#include "expand.H"
#include <Basis.H>

// Machine constant for Legendre
//
constexpr double MINEPS = 20.0*std::numeric_limits<double>::min();

Basis::Basis(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{
  use_external = false;
}

void Basis::legendre_R(int lmax, double x, Eigen::MatrixXd& p)
{
  double fact, somx2, pll, pl1, pl2;

  p(0, 0) = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (int m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p(m, m) = pll;
      if (std::isnan(p(m, m)))
	cerr << "legendre_R: p[" << m << "][" << m << "]: pll=" << pll << "\n";
      fact += 2.0;
    }
  }

  for (int m=0; m<lmax; m++) {
    pl2 = p(m, m);
    p(m+1, m) = pl1 = x*(2*m+1)*pl2;
    for (int l=m+2; l<=lmax; l++) {
      p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      if (std::isnan(p(l, m)))
	cerr << "legendre_R: p[" << l << "][" << m << "]: pll=" << pll << "\n";

      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (std::isnan(x))
    cerr << "legendre_R: x\n";
  for(int l=0; l<=lmax; l++)
    for (int m=0; m<=l; m++)
      if (std::isnan(p(l, m)))
	cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
	     << lmax << "\n";

}

void Basis::dlegendre_R(int lmax, double x,
			Eigen::MatrixXd &p, Eigen::MatrixXd &dp)
{
  double fact, somx2, pll, pl1, pl2;
  int m, l;

  p(0, 0) = pll = 1.0;
  if (lmax > 0) {
    somx2 = sqrt( (1.0 - x)*(1.0 + x) );
    fact = 1.0;
    for (m=1; m<=lmax; m++) {
      pll *= -fact*somx2;
      p(m, m) = pll;
      fact += 2.0;
    }
  }

  for (m=0; m<lmax; m++) {
    pl2 = p(m, m);
    p(m+1, m) = pl1 = x*(2*m+1)*pl2;
    for (l=m+2; l<=lmax; l++) {
      p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
      pl2 = pl1;
      pl1 = pll;
    }
  }

  if (1.0-fabs(x) < MINEPS) {
    if (x>0) x =   1.0 - MINEPS;
    else     x = -(1.0 - MINEPS);
  }

  somx2 = 1.0/(x*x - 1.0);
  dp(0, 0) = 0.0;
  for (l=1; l<=lmax; l++) {
    for (m=0; m<l; m++)
      dp(l, m) = somx2*(x*l*p(l, m) - (l+m)*p(l-1, m));
    dp(l, l) = somx2*x*l*p(l, l);
  }
}

void Basis::sinecosine_R(int mmax, double phi,
			 Eigen::VectorXd& c, Eigen::VectorXd& s)
{
  int m;

  c[0] = 1.0;
  s[0] = 0.0;

  if (mmax>0) {
    c[1] = cos(phi);
    s[1] = sin(phi);

    for (m=2; m<=mmax; m++) {
      c[m] = 2.0*c[1]*c[m-1] - c[m-2];
      s[m] = 2.0*c[1]*s[m-1] - s[m-2];
    }
  }
}


