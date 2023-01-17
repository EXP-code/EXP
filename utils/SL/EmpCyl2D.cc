#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <yaml-cpp/yaml.h>

#include <EXPException.H>
#include <EmpCyl2D.H>

std::string EmpCyl2D::cache_name_2d = ".eof_cache_2d";

// Clutton-Brock disk potential
//
double EmpCyl2D::CB_get_potl(int M, int N, double r)
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

// Clutton-Brock disk density
//
double EmpCyl2D::CB_get_dens(int M, int N, double r)
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

// Clutton-Brock disk force
//
double EmpCyl2D::CB_get_dpot(int M, int N, double r)
{
  double r2   = r*r;
  double fac  = 1.0/(1.0 + r2);
  double fac1 = (r2 - 1.0)*fac;
  double cur0 = sqrt(fac), rcum = 1.0;

  Eigen::MatrixXd p(M+2, N+1), dp(M+1, N+1);

  // Phi recursion relation
  //
  // Need up to M+1 for dPhi/dr recursion
  //
  for (int l=0; l<=M+1; l++) {
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


  // dPhi/dR recursion
  //
  for (int l=0; l<=M; l++) {

    dp(l, 0) = p(l, 0)*l/r - p(l+1, 0);
    if (N<1) break;

    dp(l, 1) = p(l, 1)*l/r - (p(l+1, 1) - 2*p(l+1, 0));
    if (N<2) break;

    for (int nn=1; nn<N; nn++) {
      dp(l, nn+1) = p(l, nn+1)*l/r - 
	( p(l+1, nn+1) - 2*p(l+1, nn) + p(l+1, nn-1) );
    }
  }

  return dp(M, N);
}


// Normalization for CB inner product
//
double EmpCyl2D::CB_norm(int n, int m)
{
  double ans = 1.0;
 
  for (int i=n+1; i<=n+2*m; i++) ans *= i;

  return pow(0.5, 2*m+1)*ans;
}


class EmpCyl2D::ExponCyl : public EmpCyl2D::ModelCyl
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

class EmpCyl2D::KuzminCyl : public EmpCyl2D::ModelCyl
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

class EmpCyl2D::MestelCyl : public EmpCyl2D::ModelCyl
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


std::shared_ptr<EmpCyl2D::ModelCyl>
EmpCyl2D::ModelCyl::createModel(const std::string type, double A)
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


double EmpCyl2D::Mapping::r_to_xi(double r)
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
    
double EmpCyl2D::Mapping::xi_to_r(double xi)
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

double EmpCyl2D::Mapping::d_xi_to_r(double xi)
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

EmpCyl2D::EmpCyl2D(int mmax, int nmax, int knots, int numr,
		   double rmin, double rmax, double A, double scale,
		   bool cmap, bool logr, const std::string type) :
  mmax(mmax), nmax(nmax), knots(knots), numr(numr),
  rmin(rmin), rmax(rmax), A(A), scale(scale), cmap(cmap), logr(logr),
  model(type)
{
  if (not read_cached_tables()) create_tables();
}


void EmpCyl2D::create_tables()
{
  Mapping  map(scale, cmap);
  auto     disk = ModelCyl::createModel(model, A);

  LegeQuad lw(knots);

  potl_array.resize(mmax+1);
  dens_array.resize(mmax+1);
  dpot_array.resize(mmax+1);
  rot_matrix.resize(mmax+1);

  Eigen::MatrixXd D(nmax, nmax);

  for (int m=0; m<=mmax; m++) {

    // Compute the covariance matrix D
    //
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
	    CB_get_potl(m, j, rr) * CB_get_potl(m, l, rr)
	    / sqrt(CB_norm(j, m)*CB_norm(l, m));
	}
      }
    }
    
    // Perform the SVD
    //
    Eigen::BDCSVD<Eigen::MatrixXd>
      svd(D, Eigen::ComputeThinU | Eigen::ComputeThinV);

    auto S = svd.singularValues();
    auto U = svd.matrixU();
    
    // Compute basis functions
    //
    double lrmin = rmin, lrmax = rmax;
    if (logr) {
      lrmin = log(rmin);
      lrmax = log(rmax);
    }
    double dr = (lrmax - lrmin)/(numr - 1), r;

    Eigen::VectorXd pot(nmax), den(nmax), dph(nmax);

    potl_array[m].resize(numr, nmax);
    dens_array[m].resize(numr, nmax);
    dpot_array[m].resize(numr, nmax);
    rot_matrix[m] = U.transpose();
    
    if (m==0) xgrid.resize(numr);

    for (int i=0; i<numr; i++) {
      r = lrmin + dr*i;
      if (logr) r = exp(r);
      if (m==0) xgrid[i] = r;

      for (int n=0; n<nmax; n++) {
	pot(n) = CB_get_potl(m, n, r) / sqrt(CB_norm(n, m));
	den(n) = CB_get_dens(m, n, r) / sqrt(CB_norm(n, m));
	dph(n) = CB_get_dpot(m, n, r) / sqrt(CB_norm(n, m));
      }
      
      pot = U.transpose() * pot;
      den = U.transpose() * den;
      dph = U.transpose() * dph;
      
      for (int n=0; n<nmax; n++) {
	potl_array[m](i, n) = pot(n);
	dens_array[m](i, n) = den(n);
	dpot_array[m](i, n) = dph(n);
      }
    }
  }

  write_cached_tables();
}

void EmpCyl2D::writeBasis(int M, const std::string& filename)
{
  std::ofstream out(filename);
  if (out) {
    
    out << std::endl << "# EOF basis grid for M=" << M << ":" << std::endl;

    for (int i=0; i<numr; i++) {
      out << std::setw(16) << xgrid[i];
      for (int n=0; n<nmax; n++) {
	out << std::setw(16) << potl_array[M](i, n)
	    << std::setw(16) << dens_array[M](i, n)
	    << std::setw(16) << dpot_array[M](i, n);
      }
      out << std::endl;
    }
  } else {
    throw std::runtime_error("EmpCyl2D::writeBasis: error opening <" + filename + ">");
  }
}

void EmpCyl2D::writeTrans(int M, const std::string& filename)
{
  std::ofstream out(filename);
  if (out) {
    
    out << std::endl << "# EOF rotation matrix for M=" << M << ":" << std::endl;
    out << rot_matrix[M] << std::endl;

  } else {
    throw std::runtime_error("EmpCyl2D::writeTrans: error opening <" + filename + ">");
  }
}


void EmpCyl2D::orthoCheck(int M, const std::string& filename)
{
  Eigen::MatrixXd orth(nmax, nmax);
  orth.setZero();

  for (int i=1; i<numr; i++) {
    for (int j=0; j<nmax; j++) {
      for (int l=0; l<nmax; l++) {
	orth(j, l) +=
	  (xgrid[i-1] * potl_array[M](i-1, j) * dens_array[M](i-1, l) +
	   xgrid[i  ] * potl_array[M](i  , j) * dens_array[M](i  , l) ) * (xgrid[i] - xgrid[i-1]) * 0.5;
      }
    }
  }
  
  std::ofstream out(filename);
  if (out) {
    out << std::endl << "# Orthogonality matrix of EOF basis for M=" << M << ":"
	<< std::endl << std::endl << orth*2.0*M_PI << std::endl;
  } else {
    throw std::runtime_error("EmpCyl2D::orthoCheck: error opening <" + filename + ">");
  }
}


bool EmpCyl2D::read_cached_tables()
{
  std::ifstream in(cache_name_2d);
  if (!in) return false;

  int MMAX, NMAX, NUMR, KNOTS;
  double RMIN, RMAX, AA, SCL;
  bool LOGR, CMAP;
  std::string MODEL;

  std::cout << "---- EmpCyl2D::read_cached_table: trying to read cached table . . ."
	    << std::endl;

  // Attempt to read magic number
  //
  unsigned int tmagic;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

  if (tmagic == hmagic) {
    
    // YAML size
    //
    unsigned ssize;
    in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

    // Make and read char buffer
    //
    auto buf = std::make_unique<char[]>(ssize+1);
    in.read(buf.get(), ssize);
    buf[ssize] = 0;		// Null terminate
    
    YAML::Node node;
    
    try {
      node = YAML::Load(buf.get());
    }
    catch (YAML::Exception& error) {
      std::ostringstream sout;
      sout << "YAML: error parsing <" << buf.get() << "> "
	   << "in " << __FILE__ << ":" << __LINE__ << std::endl
	   << "YAML error: " << error.what() << std::endl;
      throw GenericError(sout.str(), __FILE__, __LINE__, 1042, false);
    }
    
    // Get parameters
    //
    MMAX     = node["mmax"   ].as<int>();
    NMAX     = node["nmax"   ].as<int>();
    NUMR     = node["numr"   ].as<int>();
    KNOTS    = node["knots"  ].as<int>();
    LOGR     = node["logr"   ].as<bool>();
    CMAP     = node["cmap"   ].as<bool>();
    RMIN     = node["rmin"   ].as<double>();
    RMAX     = node["rmax"   ].as<double>();
    SCL      = node["scale"  ].as<double>();
    AA       = node["A"      ].as<double>();
    MODEL    = node["model"  ].as<std::string>();
  } else {
    std::cout << "---- EmpCyl2D: bad magic number in cache file" << std::endl;
    return false;
  }
    
  if (MMAX!=mmax) {
    std::cout << "---- EmpCyl2D::read_cached_table: found mmax=" << MMAX
	      << " wanted " << mmax << std::endl;
    return false;
  }

  if (NMAX!=nmax) {
    std::cout << "---- EmpCyl2D::read_cached_table: found nmax=" << NMAX
	      << " wanted " << nmax << std::endl;
    return false;
  }

  if (NUMR!=numr) {
    std::cout << "---- EmpCyl2D::read_cached_table: found numr=" << NUMR
	      << " wanted " << numr << std::endl;
    return false;
  }

  if (KNOTS!=knots) {
    std::cout << "---- EmpCyl2D::read_cached_table: found knots=" << KNOTS
	      << " wanted " << knots << std::endl;
    return false;
  }

  if (CMAP!=cmap) {
    std::cout << "---- EmpCyl2D::read_cached_table: found cmap=" << std::boolalpha << CMAP
	      << " wanted " << std::boolalpha << cmap << std::endl;
    return false;
  }

  if (LOGR!=logr) {
    std::cout << "---- EmpCyl2D::read_cached_table: found logr=" << std::boolalpha << LOGR
	      << " wanted " << std::boolalpha << logr << std::endl;
    return false;
  }

  if (RMIN!=rmin) {
    std::cout << "---- EmpCyl2D::read_cached_table: found rmin=" << RMIN
	      << " wanted " << rmin << std::endl;
    return false;
  }

  if (RMAX!=rmax) {
    std::cout << "---- EmpCyl2D::read_cached_table: found rmax=" << RMAX
	      << " wanted " << rmax << std::endl;
    return false;
  }

  if (SCL!=scale) {
    std::cout << "---- EmpCyl2D::read_cached_table: found scale=" << SCL
	      << " wanted " << scale << std::endl;
    return false;
  }

  if (AA!=A) {
    std::cout << "---- EmpCyl2D::read_cached_table: found A=" << AA
	      << " wanted " << A << std::endl;
    return false;
  }

  if (MODEL!=model) {
    std::cout << "---- EmpCyl2D::read_cached_table: found ID=" << MODEL
	      << " wanted " << model << std::endl;
    return false;
  }


  potl_array.resize(mmax+1);
  dens_array.resize(mmax+1);
  dpot_array.resize(mmax+1);
  rot_matrix.resize(mmax+1);

  xgrid.resize(numr);
  in.read((char *)xgrid.data(), numr*sizeof(double));

  for (int m=0; m<=mmax; m++) {

    potl_array[m].resize(numr, nmax);
    dens_array[m].resize(numr, nmax);
    dpot_array[m].resize(numr, nmax);
    rot_matrix[m].resize(nmax, nmax);
    
    in.read((char *)potl_array[m].data(), potl_array[m].size()*sizeof(double));
    in.read((char *)dens_array[m].data(), dens_array[m].size()*sizeof(double));
    in.read((char *)dpot_array[m].data(), dpot_array[m].size()*sizeof(double));
    in.read((char *)rot_matrix[m].data(), rot_matrix[m].size()*sizeof(double));
  }

  std::cout << "---- EmpCyl2D::read_cached_table: success!" << std::endl;

  return true;
}

void EmpCyl2D::write_cached_tables()
{
  std::ofstream out(cache_name_2d);
  if (!out) {
    std::cerr << "EmpCyl2D: error writing <" << cache_name_2d << ">" << std::endl;
    return;
  }

  // This is a node of simple {key: value} pairs.  More general
  // content can be added as needed.
  YAML::Node node;

  node["mmax"   ] = mmax;
  node["nmax"   ] = nmax;
  node["numr"   ] = numr;
  node["knots"  ] = knots;
  node["logr"   ] = logr;
  node["cmap"   ] = cmap;
  node["rmin"   ] = rmin;
  node["rmax"   ] = rmax;
  node["scale"  ] = scale;
  node["A"      ] = A;
  node["model"  ] = model;
    
  // Serialize the node
  //
  YAML::Emitter y; y << node;
  
  // Get the size of the string
  //
  unsigned int hsize = strlen(y.c_str());
  
  // Write magic #
  //
  out.write(reinterpret_cast<const char *>(&hmagic),   sizeof(unsigned int));

  // Write YAML string size
  //
  out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));
  
  // Write YAML string
  //
  out.write(reinterpret_cast<const char *>(y.c_str()), hsize);

  // Now, write the tables
  //
  out.write((char *)xgrid.data(), numr*sizeof(double));

  for (int m=0; m<=mmax; m++) {
    out.write((char *)potl_array[m].data(), potl_array[m].size()*sizeof(double));
    out.write((char *)dens_array[m].data(), dens_array[m].size()*sizeof(double));
    out.write((char *)dpot_array[m].data(), dpot_array[m].size()*sizeof(double));
    out.write((char *)rot_matrix[m].data(), rot_matrix[m].size()*sizeof(double));
  }

  std::cout << "---- EmpCyl2D::write_cached_table: cache written" << std::endl;

  return;
}
