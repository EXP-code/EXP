#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <gaussQ.H>
#include <cxxopts.H>

//! For comparing CB to 3d pillbox with Kuzmin disk and checking CB
//! orthogonality [Potential]
double CB_get_potl(int M, int N, double r)
{
  double r2   = r*r;
  double fac  = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur0 = sqrt(fac), rcum = 1.0;

  Eigen::MatrixXd p(M+1, N+1);

  for (int l=0; l<=M; l++) {
    double cur = cur0;

    p(l, 0) = cur*rcum;
    double curl1 = 0.0;

    for (int nn=0; nn<N; nn++) {
      double curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/(nn+1))*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/(nn+1))*curl2;
      p(l, nn+1) = cur*rcum;
    }
    cur0 *= fac*(2*(l+1) - 1);
    rcum *= r;
  }

  return p(M, N);
}

//! For comparing CB to 3d pillbox with Kuzmin disk and checking CB
//! orthogonality [Density]
double CB_get_dens(int M, int N, double r)
{
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur, curl1, curl2, cur0 = sqrt(fac);
  double rcum = 0.5/M_PI;
  
  Eigen::MatrixXd d(M+1, N+1), w(M+2, N+1);

  for (int l=0; l<=M+1; l++) {
    cur = cur0;

    w(l, 0) = cur;
    curl1 = 0.0;

    for (int nn=0; nn<N; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*l-1)/(nn+1))*fac1*curl1 - 
	(1.0 + (double)(2*l-1)/(nn+1))*curl2;
      w(l, nn+1) = cur;
    }
    cur0 *= fac*(2*(l+1) - 1);
  }

  for (int l=0; l<=M; l++) {
    d(l, 0) = w(l+1, 0)*rcum;
    if (N>0) {
      d(l, 1) = w(l+1, 1)*rcum;
      for (int nn=1; nn<N; nn++)
	d(l, nn+1) = (w(l+1, nn+1) - w(l+1, nn-1))*rcum;
    }

    rcum *= r;
  }

  return d(M, N);
}

//! Normalization for CB inner product
double CB_norm(int n, int m)
{
  double ans = 1.0;
 
  for (int i=n+1; i<=n+2*m; i++) ans *= i;

  return pow(0.5, 2*m+1)*ans;
}


class ModelCyl
{

protected:

  double A;
  std::string id;

public:

  static std::shared_ptr<ModelCyl>
  createModel(const std::string type, double A);

  virtual double pot (double r) = 0;
  virtual double dpot(double r) = 0;
  virtual double dens(double r) = 0;
};


class ExponCyl : public ModelCyl
{

public:

  ExponCyl(double scl) { A = scl; id = "expon"; }

  double pot(double r) {
    double y = 0.5 * r / A;
    return -2.0*M_PI*A*y*
      (std::cyl_bessel_i(0, y)*std::cyl_bessel_k(1, y) -
       std::cyl_bessel_i(1, y)*std::cyl_bessel_k(0, y));
  }

  double dpot(double r) {
    double y = 0.5 * r / A;
   return 4.0*M_PI*A*y*y*
     (std::cyl_bessel_i(0, y)*std::cyl_bessel_k(0, y) -
      std::cyl_bessel_i(1, y)*std::cyl_bessel_k(1, y));
  }

  double dens(double r) {
    // This 4pi from Poisson's eqn
    //        |
    //        |       /-- This begins the true projected density profile
    //        |       |
    //        v       v
    return 4.0*M_PI * exp(-r/A);
  }

};

class KuzminCyl : public ModelCyl
{
public:
  
  KuzminCyl(double scl) { A = scl; id = "kuzmin"; }

  double pot(double R) {
    double a2 = A * A;
    return -1.0/sqrt(R*R + a2);
  }

  double dpot(double R) {
    double a2 = A * A;
    return R/pow(R*R + a2, 1.5);
  }

  double dens(double R) {
    double a2 = A * A;
    return 4.0*M_PI*A/pow(R*R + a2, 1.5)/(2.0*M_PI);
    //     ^
    //     |
    // This 4pi from Poisson's eqn
  }

};

class MestelCyl : public ModelCyl
{
public:
  
  MestelCyl(double scl) { A = scl;  id = "mestel"; }

  double pot(double R) {
    return M_PI/(2.0*A)*log(0.5*R/A);
  }

  double dpot(double R) {
    double a2 = A * A;
    double fac = sqrt(1.0 + R*R/a2);
    return M_PI/(2.0*A*R);
  }

  double dens(double R) {
    if (R>A)
      return 0.0;
    else
      return 4.0*M_PI/(2.0*M_PI*A*R)*acos(R/A);
      //     ^
      //     |
      // This 4pi from Poisson's eqn
  }
};


std::shared_ptr<ModelCyl>
ModelCyl::createModel(const std::string type, double A)
{
  // Convert ID string to lower case
  //
  std::string data(type);
  std::transform(data.begin(), data.end(), data.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  if (data.find("kuzmin") != std::string::npos) {
    std::cout << "Making a Kuzmin disk" << std::endl;
    return std::make_shared<KuzminCyl>(A);
  }

  if (data.find("mestel") != std::string::npos) {
    std::cout << "Making a finite Mestel disk" << std::endl;
    return std::make_shared<MestelCyl>(A);
  }

  if (data.find("expon") != std::string::npos) {
    std::cout << "Making an Exponential disk" << std::endl;
    return std::make_shared<ExponCyl>(A);
  }

  // Default if nothing else matches
  std::cout << "[Default] making an Exponential disk" << std::endl;
  return std::make_shared<ExponCyl>(A);
}


class Mapping
{
protected:
  bool cmap;
  double scale;

public:

  Mapping(double scale, bool cmap) : scale(scale), cmap(cmap) {}

  double r_to_xi(double r)
  {
    if (r<0.0) {
      std::ostringstream ostr;
      ostr << "radius=" << r << " < 0!";
      throw std::runtime_error(ostr.str());
    }

    if (cmap) {
      return (r/scale-1.0)/(r/scale+1.0);
    } else {
      return r;
    }
  }
    
  double xi_to_r(double xi)
  {
    if (cmap) {
      if (xi<-1.0) {
	std::ostringstream ostr;
	ostr << "xi=" << xi << " < -1!";
	throw std::runtime_error(ostr.str());
      }

      if (xi>=1.0) {
	std::ostringstream ostr;
	ostr << "xi=" << xi << " >= 1!";
	throw std::runtime_error(ostr.str());
      }

      return (1.0+xi)/(1.0 - xi) * scale;
    } else {
      return xi;
    }
  }

  double d_xi_to_r(double xi)
  {
    if (cmap) {
      if (xi<-1.0) {
	std::ostringstream ostr;
	ostr << "xi=" << xi << " < -1!";
	throw std::runtime_error(ostr.str());
      }

      if (xi>=1.0) {
	std::ostringstream ostr;
	ostr << "xi=" << xi << " >= 1!";
	throw std::runtime_error(ostr.str());
      }
    
      return 0.5*(1.0-xi)*(1.0-xi)/scale;
    } else {
      return 1.0;
    }
  }
};


int main(int argc, char** argv)
{
  bool logr = false, cmap = false, ortho = false, trans = false;
  int numr, mmax, nmax, knots, num, M;
  double A, scale, rmin, rmax;
  std::string filename, type;

  // Parse command line
  //
  cxxopts::Options options(argv[0], "Check the consistency a spherical SL basis");

  options.add_options()
    ("h,help", "Print this help message")
    ("logr", "Plot output grid with logarithmic spacing")
    ("cmap", "Use mapped coordinates")
    ("ortho", "Compute EOF orthogonal matrix")
    ("trans", "Print the rotation matrix")
    ("scale", "scaling from real coordinates to table",
     cxxopts::value<double>(scale)->default_value("1.0"))
    ("M,harmonic", "Aximuthal harmonic",
     cxxopts::value<int>(M)->default_value("0"))
    ("A,length", "characteristic disk scale length",
     cxxopts::value<double>(A)->default_value("1.0"))
    ("mmax", "maximum number of angular harmonics in the expansion",
     cxxopts::value<int>(mmax)->default_value("4"))
    ("nmax", "maximum number of radial harmonics in the expansion",
     cxxopts::value<int>(nmax)->default_value("10"))
    ("numr", "radial knots for the SL grid",
     cxxopts::value<int>(numr)->default_value("1000"))
    ("r,rmin", "minimum radius for the SL grid",
     cxxopts::value<double>(rmin)->default_value("0.0001"))
    ("R,rmax", "maximum radius for the SL grid",
     cxxopts::value<double>(rmax)->default_value("20.0"))
    ("num", "number of output grid points",
     cxxopts::value<int>(num)->default_value("1000"))
    ("knots", "Number of Legendre integration knots",
     cxxopts::value<int>(knots)->default_value("40"))
    ("type", "Target model type (kuzmin, mestel, expon)",
     cxxopts::value<std::string>(type)->default_value("expon"))
    ("o,filename", "Output filename",
     cxxopts::value<std::string>(filename)->default_value("testeof.dat"))
    ;

  
  //===================
  // Parse options
  //===================

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return 2;
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    std::cout << options.help() << std::endl << std::endl;
    return 1;
  }

  // Log spacing?
  //
  if (vm.count("logr")) logr = true;

  // Mapped coordinates?
  //
  if (vm.count("cmap")) cmap = true;

  // Print rotation matrix?
  //
  if (vm.count("trans")) trans = true;

  // Check orthogonality?
  //
  if (vm.count("ortho")) ortho = true;

  Mapping  map(1.0, true);
  auto     disk = ModelCyl::createModel(type, A);

  std::ofstream out(filename);
  if (!out) {
    std::cout << "Can't open <" << filename << "> for output" << endl;
    exit(-1);
  }

  LegeQuad lw(knots);

  // Compute the covariance matrix D
  //
  Eigen::MatrixXd D(nmax, nmax);
  D.setZero();

  double ximin = map.r_to_xi(rmin);
  double ximax = map.r_to_xi(rmax);

  for (int k=0; k<knots; k++) {
    double xx  = ximin + (ximax - ximin)*lw.knot(k);
    double rr  = map.xi_to_r(xx);
    double fac = lw.weight(k) * rr / map.d_xi_to_r(xx) * (ximax - ximin);
    
    for (int j=0; j<nmax; j++) {
      for (int l=0; l<nmax; l++) {
	D(j, l) += fac * disk->dens(rr) *
	  CB_get_potl(M, j, rr) * CB_get_potl(M, l, rr)
	  / sqrt(CB_norm(j, M)*CB_norm(l, M));
      }
    }
  }

  // Perform the SVD
  //
  
  Eigen::BDCSVD<Eigen::MatrixXd>
    svd(D, Eigen::ComputeThinU | Eigen::ComputeThinV);
  auto S = svd.singularValues();
  auto U = svd.matrixU();

  std::cout << std::endl << "Eigenvalues: " << std::endl << S << std::endl;

  if (trans) {
    std::cout << std::endl << "Transformation matrix: " << std::endl
	      << std::endl << U.transpose() << std::endl;
  }

  // Compute some basis functions
  //
  double lrmin = rmin, lrmax = rmax;
  if (logr) {
    lrmin = log(rmin);
    lrmax = log(rmax);
  }
  double dr = (lrmax - lrmin)/(numr - 1), r;

  Eigen::VectorXd pot(nmax), den(nmax);

  for (int i=0; i<numr; i++) {
    r = lrmin + dr*i;
    if (logr) r = exp(r);

    for (int n=0; n<nmax; n++) {
      pot(n) = CB_get_potl(M, n, r) / sqrt(CB_norm(n, M));
      den(n) = CB_get_dens(M, n, r) / sqrt(CB_norm(n, M));
    }

    pot = U.transpose() * pot;
    den = U.transpose() * den;

    out << std::setw(16) << r;
    for (int n=0; n<nmax; n++) {
      out << std::setw(16) << pot(n)
	  << std::setw(16) << den(n);
    }
    out << std::endl;
  }

  if (ortho) {

    Eigen::MatrixXd orth(nmax, nmax);
    orth.setZero();

    for (int k=0; k<knots; k++) {
      double xx  = ximin + (ximax - ximin)*lw.knot(k);
      double rr  = map.xi_to_r(xx);
      double fac = lw.weight(k) * rr / map.d_xi_to_r(xx) * (ximax - ximin);
      
      for (int n=0; n<nmax; n++) {
	pot(n) = CB_get_potl(M, n, rr) / sqrt(CB_norm(n, M));
	den(n) = CB_get_dens(M, n, rr) / sqrt(CB_norm(n, M));
      }
      
      pot = U.transpose() * pot;
      den = U.transpose() * den;
      
      for (int j=0; j<nmax; j++) {
	for (int l=0; l<nmax; l++) {
	  orth(j, l) += fac * pot(j) * den(l);
	}
      }
    }
    
    std::cout << std::endl << "Orthogonality matrix of EOF basis:"
	      << std::endl << std::endl << orth*2.0*M_PI << std::endl;
  }

  return 0;
}
