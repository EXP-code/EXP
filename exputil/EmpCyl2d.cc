#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <cmath>

#include <yaml-cpp/yaml.h>

#include <EXPException.H>
#include <localmpi.H>
#include <EmpCyl2d.H>
#include <EXPmath.H>

// Set to true for orthogonality checking
//
bool EmpCyl2d::Basis2d::debug = false;

// Clutton-Brock two-dimensional disk
//
class EmpCyl2d::CluttonBrock : public EmpCyl2d::Basis2d
{
protected:

  void test_ortho();

public:

  CluttonBrock() { test_ortho(); }

  double potl(int M, int N, double r);
  double dens(int M, int N, double r);
  double dpot(int M, int N, double r);
  double norm(int N, int M);
};

Eigen::VectorXd bessjz(int n, int m);

// Clutton-Brock two-dimensional disk
//
class EmpCyl2d::Bessel : public EmpCyl2d::Basis2d
{
private:
  std::vector<Eigen::VectorXd> roots, normv;
  int mmax, nmax;
  double L;

  void test_ortho();

public:

  Bessel(int mmax, int nmax, double L);

  double potl(int M, int N, double r);
  double dens(int M, int N, double r);
  double dpot(int M, int N, double r);
  double norm(int N, int M);
};

EmpCyl2d::Bessel::Bessel(int mmax, int nmax, double L) :
  mmax(mmax), nmax(nmax), L(L)
{
  roots.resize(mmax+1);
  normv.resize(mmax+1);
  for (int m=0; m<=mmax; m++) {
    roots[m] = bessjz(m, nmax);
    normv[m].resize(nmax);
    for (int n=0; n<nmax; n++) {
      double eval = EXPmath::cyl_bessel_j(m+1, roots[m][n]) * L;
      normv[m][n] = eval*eval*0.5 * roots[m][n]/L;
    }
  }

  test_ortho();
}

void EmpCyl2d::Bessel::test_ortho()
{
  // Sanity check; set to false for production
  //
  if (debug and myid==0) {

    const int knots = 128;
    LegeQuad lw(knots);
    Eigen::MatrixXd orth(nmax, nmax);

    std::cout << std::endl;
    for (int m=0; m<=mmax; m++) {
      orth.setZero();
      for (int i=0; i<knots; i++) {
	double r = L*lw.knot(i);
	for (int j=0; j<nmax; j++) {
	  for (int l=0; l<nmax; l++) {
	    orth(j, l) += lw.weight(i)*L * r *
	      potl(m, j, r)*dens(m, l, r) / sqrt(norm(j, m)*norm(l, m));
	  }
	}
      }
      std::cout << "# Bessel orthgonality for M=" << m << std::endl;
      std::cout << orth*2.0*M_PI << std::endl << std::endl;
    }
  }
}

void EmpCyl2d::CluttonBrock::test_ortho()
{
  // Sanity check; set to false for production
  //
  if (debug and myid==0) {

    const int knots = 128, nmax = 10, mmax = 4;
    LegeQuad lw(knots);
    Eigen::MatrixXd orth(nmax, nmax);

    auto xi_to_r   = [](double xi) { return (1.0+xi)/(1.0 - xi);   };
    auto d_xi_to_r = [](double xi) { return 0.5*(1.0-xi)*(1.0-xi); };

    for (int m=0; m<=mmax; m++) {
      orth.setZero();
      for (int i=0; i<knots; i++) {
	double xi = 2.0*(lw.knot(i) - 0.5);
	double r = xi_to_r(xi);
	for (int j=0; j<nmax; j++) {
	  for (int l=0; l<nmax; l++) {
	    orth(j, l) += lw.weight(i)*2.0 * r / d_xi_to_r(xi) *
	      potl(m, j, r)*dens(m, l, r) / sqrt(norm(j, m)*norm(l, m));
	  }
	}
      }
      std::cout << "# CluttonBrock orthgonality for M=" << m << std::endl;
      std::cout << orth*2.0*M_PI << std::endl << std::endl;
    }
  }
}

double EmpCyl2d::CluttonBrock::potl(int M, int N, double r)
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
double EmpCyl2d::CluttonBrock::dens(int M, int N, double r)
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
double EmpCyl2d::CluttonBrock::dpot(int M, int N, double r)
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
double EmpCyl2d::CluttonBrock::norm(int n, int m)
{
  double ans = 1.0;
 
  for (int i=n+1; i<=n+2*m; i++) ans *= i;

  return pow(0.5, 2*m+1)*ans;
}


// Bessel disk potential
double EmpCyl2d::Bessel::potl(int M, int N, double r)
{
  return EXPmath::cyl_bessel_j(M, r*roots[M][N]/L);
}

// Bessel disk density
//
double EmpCyl2d::Bessel::dens(int M, int N, double r)
{
  return EXPmath::cyl_bessel_j(M, r*roots[M][N]/L)*roots[M][N]/L/(2.0*M_PI);
}

// Bessel disk force
//
double EmpCyl2d::Bessel::dpot(int M, int N, double r)
{
  double z = r*roots[M][N]/L;
  if (M==0) {		    // Using the derivative recursion relation
    return EXPmath::cyl_bessel_j(M+1, z) * roots[M][N]/L;
  } else {
    return -(EXPmath::cyl_bessel_j(M-1, z) - EXPmath::cyl_bessel_j(M+1, z))
      * 0.5*roots[M][N]/L;
  }
}

// Normalization for Bessel inner product
//
double EmpCyl2d::Bessel::norm(int n, int m)
{
  return normv[m][n];
}

std::shared_ptr<EmpCyl2d::Basis2d>
EmpCyl2d::Basis2d::createBasis(int mmax, int nmax, double rmax,
			       const std::string& type)
{
  // Convert ID string to lower case
  //
  std::string data(type);
  std::transform(data.begin(), data.end(), data.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  if (data.find("cb") != std::string::npos) {
    if (myid==0)
      std::cout << "---- EmpCyl2d::Basis2d: Making a Clutton-Brock basis" << std::endl;
    return std::make_shared<CluttonBrock>();
  }

  if (data.find("bess") != std::string::npos) {
    if (myid==0)
      std::cout << "---- EmpCyl2d::Basis2d: Making a finite Bessel basis" << std::endl;
    return std::make_shared<Bessel>(mmax, nmax, rmax);
  }

  // Default if nothing else matches
  if (myid==0)
    std::cout << "---- EmpCyl2d::Basis2d: Making an Clutton-Brock basis [Default]" << std::endl;
  return std::make_shared<CluttonBrock>();
}

class EmpCyl2d::ExponCyl : public EmpCyl2d::ModelCyl
{

private:

  double acyl, sigma0;

  // Assign values from YAML
  //
  void parse(const YAML::Node& conf)
  {
    try {
      if (conf["acyl"]) 
	acyl = conf["acyl"].as<double>();
      else
	acyl = 1.0;
    }
    catch (YAML::Exception & error) {
      if (myid==0)
	std::cout << "Error parsing parameters in EmpCyl2d::ExponCyl: "
		  << error.what() << std::endl
		  << std::string(60, '-') << std::endl
		  << "Config node"        << std::endl
		  << std::string(60, '-') << std::endl
		  << conf                 << std::endl
		  << std::string(60, '-') << std::endl;
      MPI_Finalize();
      exit(-1);
    }
  }

public:

  ExponCyl(const YAML::Node& par)
  {
    parse(par);
    sigma0 = 0.5/(M_PI*acyl*acyl);
    id = "expon";
  }

  double pot(double r) {
    double y = 0.5 * r / acyl;
    return -M_PI*sigma0*r *
      (EXPmath::cyl_bessel_i(0, y)*EXPmath::cyl_bessel_k(1, y) -
       EXPmath::cyl_bessel_i(1, y)*EXPmath::cyl_bessel_k(0, y));
  }

  double dpot(double r) {
    double y = 0.5 * r / acyl;
    return 2.0*M_PI*sigma0*y*
     (EXPmath::cyl_bessel_i(0, y)*EXPmath::cyl_bessel_k(0, y) -
      EXPmath::cyl_bessel_i(1, y)*EXPmath::cyl_bessel_k(1, y));
  }

  double dens(double r) {
    return sigma0*exp(-r/acyl);
  }

};

class EmpCyl2d::KuzminCyl : public EmpCyl2d::ModelCyl
{
private:

  double acyl;

  // Assign values from YAML
  //
  void parse(const YAML::Node& conf)
  {
    try {
      if (conf["acyl"]) 
	acyl = conf["acyl"].as<double>();
      else
	acyl = 1.0;
    }
    catch (YAML::Exception & error) {
      if (myid==0)
	std::cout << "Error parsing parameters in EmpCyl2d::KuzminCyl: "
		  << error.what() << std::endl
		  << std::string(60, '-') << std::endl
		  << "Config node"        << std::endl
		  << std::string(60, '-') << std::endl
		  << conf                 << std::endl
		  << std::string(60, '-') << std::endl;
      MPI_Finalize();
      exit(-1);
    }
  }

public:
  
  KuzminCyl(const YAML::Node& par)
  {
    parse(par);
    id = "kuzmin";
  }

  double pot(double R) {
    double a2 = acyl * acyl;
    return -1.0/sqrt(R*R + a2);
  }

  double dpot(double R) {
    double a2 = acyl * acyl;
    return R/pow(R*R + a2, 1.5);
  }

  double dens(double R) {
    double a2 = acyl * acyl;
    return 4.0*M_PI*acyl/pow(R*R + a2, 1.5)/(2.0*M_PI);
    //     ^
    //     |
    // This 4pi from Poisson's eqn
  }

};


class EmpCyl2d::MestelCyl : public EmpCyl2d::ModelCyl
{
protected:

  double vrot, rot;

  // Assign values from YAML
  //
  void parse(const YAML::Node& conf)
  {
    try {
      if (conf["vrot"]) 
	vrot = conf["vrot"].as<double>();
      else
	vrot = 1.0;
    }
    catch (YAML::Exception & error) {
      if (myid==0)
	std::cout << "Error parsing parameters in EmpCyl2d::MestelCyl: "
		  << error.what() << std::endl
		  << std::string(60, '-') << std::endl
		  << "Config node"        << std::endl
		  << std::string(60, '-') << std::endl
		  << conf                 << std::endl
		  << std::string(60, '-') << std::endl;
      MPI_Finalize();
      exit(-1);
    }
  }

public:
  
  MestelCyl(const YAML::Node& par)
  {
    parse(par);
    rot = vrot*vrot;
    id = "mestel";
  }

  virtual double pot(double R) {
    if (R>0.0) return rot*log(R);
    else throw std::runtime_error("MestelCyl::pot: R<=0");
  }

  virtual double dpot(double R) {
    if (R>0.0) return rot/R;
    else throw std::runtime_error("MestelCyl::dpot: R<=0");
  }

  virtual double dens(double R) {
    if (R>0.0) return rot/(2.0*M_PI*R);
    else throw std::runtime_error("MestelCyl::dens: R<=0");
  }
};


class EmpCyl2d::ZangCyl : public EmpCyl2d::MestelCyl
{
  
private:
  //! Parameters
  double vr, mu, nu, ri, ro;

  //! Softening factor
  double asoft = 1.0e-8;

  //! Ignore inner cut-off for N<0.05
  bool Inner = true;

  //! Taper factors
  double Tifac, Tofac;

  //! Inner taper function
  double Tinner(double Jp)
  {
    double fac = pow(Jp, nu);
    return fac/(Tifac + fac);
  }

  //! Outer taper function
  double Touter(double Jp)
  {
    return 1.0/(1.0 + pow(Jp/vr, mu));
  }

  //! Deriv of inner taper function
  double dTinner(double Jp)
  {
    double fac  = pow(Jp, nu);
    double fac2 = Tifac + fac;
    return Tifac*nu/Jp/fac2;
  }

  //! Deriv of outer taper function
  double dTouter(double Jp)
  {
    double fac = pow(Jp/Tofac, mu);
    double fac2 = 1.0 + fac;
    return -nu*fac/Jp/fac2;
  }

protected:

  //! Assign values from YAML
  void parse(const YAML::Node& conf)
  {
    try {
      if (conf["vrot"]) 
	vrot = conf["vrot"].as<double>();
      else
	vrot = 1.0;

      if (conf["Ninner"]) 
	nu = conf["Ninner"].as<double>();
      else
	nu = 2.0;

      if (conf["Mouter"]) 
	mu = conf["Mouter"].as<double>();
      else
	mu = 2.0;

      if (conf["Ri"]) 
	ri = conf["Ri"].as<double>();
      else
	ri = 1.0;

      if (conf["Ro"]) 
	ro = conf["Ro"].as<double>();
      else
	ro = 10.0;
    }
    catch (YAML::Exception & error) {
      if (myid==0)
	std::cout << "Error parsing parameters in EmpCyl2d::ZangCyl: "
		  << error.what() << std::endl
		  << std::string(60, '-') << std::endl
		  << "Config node"        << std::endl
		  << std::string(60, '-') << std::endl
		  << conf                 << std::endl
		  << std::string(60, '-') << std::endl;
      MPI_Finalize();
      exit(-1);
    }
  }

public:
  
  //! Constructor
  ZangCyl(const YAML::Node& par) : MestelCyl(par)
  {
    // Parse the YAML
    parse(par);
    // Assign the id
    id = "zang";

    // Cache taper factors
    Tifac = pow(ri*vr, nu);
    Tofac = ro*vr;

    if (nu<0.05) {
      // Exponent is now for mapping only
      Inner = false;
    }
  }

  //! Surface density
  double dens(double R)
  {
    double ret = MestelCyl::dens(R) * Touter(R);
    if (Inner) ret *= Tinner(R);
    return ret;
  }

};

std::shared_ptr<EmpCyl2d::ModelCyl> EmpCyl2d::createModel()
{
  try {
    // Get model name
    //
    std::string name;
    if (Params["name"])
      name = Params["name"].as<std::string>();
    else
      throw std::runtime_error("EmpCyl2d::createModel: the 'diskconf' stanza must specify the model 'name'");

    // Convert ID string to lower case
    //
    std::transform(name.begin(), name.end(), name.begin(),
		   [](unsigned char c){ return std::tolower(c); });

  if (name.find("kuzmin") != std::string::npos) {
    if (myid==0) {
      std::cout << "---- EmpCyl2d::ModelCyl: Making a Kuzmin disk";
      if (Params["parameters"]) std::cout << " with " << Params["parameters"];
      std::cout << std::endl;
    }
    return std::make_shared<KuzminCyl>(Params["parameters"]);
  }

  if (name.find("mestel") != std::string::npos) {
    if (myid==0) {
      std::cout << "---- EmpCyl2d::ModelCyl: Making a finite Mestel disk";
      if (Params["parameters"]) std::cout << " with " << Params["parameters"];
      std::cout << std::endl;
    }
    return std::make_shared<MestelCyl>(Params["parameters"]);
  }

  if (name.find("zang") != std::string::npos) {
    if (myid==0) {
      std::cout << "---- EmpCyl2d::ModelCyl: Making a double-tapered Zang";
      if (Params["parameters"]) std::cout << " with " << Params["parameters"];
      std::cout << std::endl;
    }
    return std::make_shared<ZangCyl>(Params["parameters"]);
  }

  if (name.find("expon") != std::string::npos) {
    if (myid==0) {
      std::cout << "---- EmpCyl2d::ModelCyl: Making an Exponential disk";
      if (Params["parameters"]) std::cout << " with " << Params["parameters"];
      std::cout << std::endl;
    }
    return std::make_shared<ExponCyl>(Params["parameters"]);
  }

  // Default if nothing else matches
  if (myid==0)
    std::cout << "---- EmpCyl2d::ModelCyl: Making an Exponential disk [Default]"
	      << std::endl;
  return std::make_shared<ExponCyl>(Params["parameters"]);
  }
  catch (YAML::Exception & error) {
    if (myid==0)
      std::cout << "Error parsing parameters in EmpCyl2d::modelCreate: "
		<< error.what() << std::endl
		<< std::string(60, '-') << std::endl
		<< "Config node"        << std::endl
		<< std::string(60, '-') << std::endl
		<< Params               << std::endl
		<< std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}


double EmpCyl2d::Mapping::r_to_xi(double r)
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
    
double EmpCyl2d::Mapping::xi_to_r(double xi)
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

double EmpCyl2d::Mapping::d_xi_to_r(double xi)
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

// Default cache name
//
const std::string EmpCyl2d::default_cache_name = ".eof_cache_2d";


// The main constructor
//
EmpCyl2d::EmpCyl2d
(int mmax, int nmaxfid, int nmax, int knots, int numr,
 double rmin, double rmax, double scale, bool cmap, bool logr,
 const YAML::Node& par,
 const std::string biorth, const std::string cache) :
  mmax(mmax), nmaxfid(nmaxfid), nmax(nmax), knots(knots), numr(numr),
  rmin(rmin), rmax(rmax), scale(scale), cmap(cmap), logr(logr),
  Params(par), biorth(biorth), cache_name_2d(cache)
{
  if (cache_name_2d.size()==0) cache_name_2d = default_cache_name;

  disk  = createModel();
  basis = Basis2d::createBasis(mmax, nmaxfid, rmax, biorth);

  basis_test = false;

  if (not ReadH5Cache()) create_tables();

  configured = true;

  nmax = std::min<int>(nmax, nmaxfid);
}


EmpCyl2d::EmpCyl2d
(int mmax, int nmaxfid, int nmax, int knots, int numr,
 double rmin, double rmax, double scale, bool cmap, bool logr,
 std::shared_ptr<EmpCyl2d::ModelCyl> disk,
 const std::string biorth, const std::string cache) :
  mmax(mmax), nmaxfid(nmaxfid), nmax(nmax), knots(knots), numr(numr),
  rmin(rmin), rmax(rmax), scale(scale), cmap(cmap), logr(logr),
  disk(disk), biorth(biorth), cache_name_2d(cache)
{
  if (cache_name_2d.size()==0) cache_name_2d = default_cache_name;

  model = disk->ID();
  basis = Basis2d::createBasis(mmax, nmaxfid, rmax, biorth);

  basis_test = false;

  if (not ReadH5Cache()) create_tables();

  configured = true;

  nmax = std::min<int>(nmax, nmaxfid);
}


void EmpCyl2d::create_tables()
{
  map = Mapping(scale, cmap);

  LegeQuad lw(knots);

  potl_array.resize(mmax+1);
  dens_array.resize(mmax+1);
  dpot_array.resize(mmax+1);
  rot_matrix.resize(mmax+1);

  Eigen::MatrixXd D(nmaxfid, nmaxfid);

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
    
      for (int j=0; j<nmaxfid; j++) {
	for (int l=0; l<nmaxfid; l++) {
	  D(j, l) += fac * disk->dens(rr) *
	    basis->potl(m, j, rr) * basis->potl(m, l, rr)
	    / sqrt(basis->norm(j, m)*basis->norm(l, m));
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

    Eigen::VectorXd pot(nmaxfid), den(nmaxfid), dph(nmaxfid);

    potl_array[m].resize(numr, nmax);
    dens_array[m].resize(numr, nmax);
    dpot_array[m].resize(numr, nmax);
    rot_matrix[m] = U.transpose();
    
    if (m==0) xgrid.resize(numr);

    for (int i=0; i<numr; i++) {
      r = lrmin + dr*i;
      if (logr) r = exp(r);
      if (m==0) xgrid[i] = r;

      for (int n=0; n<nmaxfid; n++) {
	pot(n) = basis->potl(m, n, r) / sqrt(basis->norm(n, m));
	den(n) = basis->dens(m, n, r) / sqrt(basis->norm(n, m));
	dph(n) = basis->dpot(m, n, r) / sqrt(basis->norm(n, m));
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

  if (myid==0) WriteH5Cache();
}

void EmpCyl2d::writeBasis(int M, const std::string& filename)
{
  std::ofstream out(filename);

  if (out) {
    out << std::endl << "# EOF basis grid for M=" << M << ":" << std::endl;

    for (int i=0; i<nmax; i++) {
      out << std::setw(16) << xgrid[i];
      for (int n=0; n<nmax; n++) {

	out << std::setw(16) << potl_array[M](i, n)
	    << std::setw(16) << dens_array[M](i, n)
	    << std::setw(16) << dpot_array[M](i, n);
      }
      out << std::endl;
    }
  } else {
    throw std::runtime_error("EmpCyl2d::writeBasis: error opening <" + filename + ">");
  }
}

void EmpCyl2d::writeTrans(int M, const std::string& filename)
{
  std::ofstream out(filename);
  if (out) {
    
    out << std::endl << "# EOF rotation matrix for M=" << M << ":" << std::endl;
    out << rot_matrix[M] << std::endl;

  } else {
    throw std::runtime_error("EmpCyl2d::writeTrans: error opening <" + filename + ">");
  }
}


void EmpCyl2d::orthoCheck(int M, const std::string& filename)
{
  Eigen::MatrixXd orth0(nmax, nmax), orth1(nmax, nmax);
  orth0.setZero();
  orth1.setZero();

  for (int i=1; i<numr; i++) {
    for (int j=0; j<nmax; j++) {
      for (int l=0; l<nmax; l++) {
	orth0(j, l) +=
	  (xgrid[i-1] * basis->potl(M, j, xgrid[i-1]) * basis->dens(M, l, xgrid[i-1]) + 
	   xgrid[i  ] * basis->potl(M, j, xgrid[i  ]) * basis->dens(M, l, xgrid[i  ]) ) *
	  (xgrid[i] - xgrid[i-1]) * 0.5 / sqrt(basis->norm(j, M)*basis->norm(l, M));

	orth1(j, l) +=
	  (xgrid[i-1] * potl_array[M](i-1, j) * dens_array[M](i-1, l) +
	   xgrid[i  ] * potl_array[M](i  , j) * dens_array[M](i  , l) ) *
	  (xgrid[i] - xgrid[i-1]) * 0.5;
      }
    }
  }
  
  std::ofstream out(filename);
  if (out) {
    out << std::endl << "# Orthogonality matrix of basis for M=" << M << ":"
	<< std::endl << std::endl << orth0*2.0*M_PI << std::endl << std::endl
	<< std::endl << "# Orthogonality matrix of EOF basis for M=" << M << ":"
	<< std::endl << std::endl << orth1*2.0*M_PI << std::endl << std::endl;
  } else {
    throw std::runtime_error("EmpCyl2d::orthoCheck: error opening <" + filename + ">");
  }
}


std::vector<Eigen::MatrixXd> EmpCyl2d::orthoCheck()
{
  std::vector<Eigen::MatrixXd> ret(mmax+1);

  for (int M=0; M<=mmax; M++) {

    ret[M].resize(nmax, nmax);
    ret[M].setZero();

    for (int i=1; i<numr; i++) {
      for (int j=0; j<nmax; j++) {
	for (int l=0; l<nmax; l++) {
	  ret[M](j, l) +=
	    (xgrid[i-1] * potl_array[M](i-1, j) * dens_array[M](i-1, l) +
	     xgrid[i  ] * potl_array[M](i  , j) * dens_array[M](i  , l) ) *
	    (xgrid[i] - xgrid[i-1]) * 0.5;
	}
      }
    }

    ret[M] *= 2.0*M_PI;
  }

  return ret;
}


void EmpCyl2d::WriteH5Cache()
{
  if (myid) return;

  try {
    // Create a new hdf5 file or overwrite an existing file
    //
    HighFive::File file(cache_name_2d, HighFive::File::Overwrite);
    
    // Workaround for lack of HighFive boolean support
    int ilogr = 0, icmap = 0;
    if (logr) ilogr = 1;
    if (cmap) icmap = 1;

    // Serialize the config and make a string
    YAML::Emitter y; y << Params;
    std::string params(y.c_str());

    // Parameters
    //
    file.createAttribute<int>        ("mmax",   HighFive::DataSpace::From(mmax)).  write(mmax);
    file.createAttribute<int>        ("nmaxfid",   HighFive::DataSpace::From(nmaxfid)).  write(nmaxfid);
    file.createAttribute<int>        ("nmax", HighFive::DataSpace::From(nmax)).write(nmax);
    file.createAttribute<int>        ("numr",   HighFive::DataSpace::From(numr)).  write(numr);
    file.createAttribute<int>        ("knots",  HighFive::DataSpace::From(knots)). write(knots);
    file.createAttribute<int>        ("ilogr",  HighFive::DataSpace::From(ilogr)). write(ilogr);
    file.createAttribute<int>        ("icmap",  HighFive::DataSpace::From(icmap)). write(icmap);
    file.createAttribute<double>     ("rmin",   HighFive::DataSpace::From(rmin)).  write(rmin);
    file.createAttribute<double>     ("rmax",   HighFive::DataSpace::From(rmax)).  write(rmax);
    file.createAttribute<double>     ("scale",  HighFive::DataSpace::From(scale)). write(scale);
    file.createAttribute<std::string>("params", HighFive::DataSpace::From(params)).write(params);
    file.createAttribute<std::string>("model",  HighFive::DataSpace::From(model)). write(model);
    file.createAttribute<std::string>("biorth", HighFive::DataSpace::From(biorth)).write(biorth);
      
    // Arrays
    //
    file.createDataSet("xgrid", xgrid);

    for (int m=0; m<=mmax; m++) {
      std::ostringstream sout;
      sout << m;
      auto harmonic = file.createGroup(sout.str());
      
      harmonic.createDataSet("potl", potl_array[m]);
      harmonic.createDataSet("dens", dens_array[m]);
      harmonic.createDataSet("dpot", dpot_array[m]);
      harmonic.createDataSet("rot",  rot_matrix[m]);
    }

  } catch (HighFive::Exception& err) {
    std::cerr << err.what() << std::endl;
  }
    
  std::cout << "---- EmpCyl2d::WriteH5Cache: "
	    << "wrote <" << cache_name_2d << ">" << std::endl;
}

bool EmpCyl2d::ReadH5Cache()
{
  try {
    // Silence the HDF5 error stack
    //
    HighFive::SilenceHDF5 quiet;
    

    // Open the hdf5 file
    //
    HighFive::File file(cache_name_2d, HighFive::File::ReadOnly);
    

    // Try checking the rest of the parameters before reading arrays
    //
    auto checkInt = [&file](int value, std::string name)
    {
      int v; HighFive::Attribute vv = file.getAttribute(name); vv.read(v);
      if (value == v) return true; return false;
    };

    auto checkDbl = [&file](double value, std::string name)
    {
      double v; HighFive::Attribute vv = file.getAttribute(name); vv.read(v);
      if (fabs(value - v) < 1.0e-16) return true; return false;
    };

    auto checkStr = [&file](std::string value, std::string name)
    {
      std::string v; HighFive::Attribute vv = file.getAttribute(name); vv.read(v);
      if (value.compare(v)==0) return true; return false;
    };

    //

    // Serialize the config and make a string for checking
    YAML::Emitter y; y << Params;
    std::string params(y.c_str());

    // Workaround for lack of HighFive boolean support
    int ilogr = 0, icmap = 0;
    if (logr) ilogr = 1;
    if (cmap) icmap = 1;
    
    if (not checkInt(mmax,     "mmax"))      return false;
    if (not checkInt(nmaxfid,  "nmaxfid"))   return false;
    if (not checkInt(nmax,     "nmax"))      return false;
    if (not checkInt(numr,     "numr"))      return false;
    if (not checkInt(knots,    "knots"))     return false;
    if (not checkInt(ilogr,    "ilogr"))     return false;
    if (not checkInt(icmap,    "icmap"))     return false;
    if (not checkDbl(rmin,     "rmin"))      return false;
    if (not checkDbl(rmax,     "rmax"))      return false;
    if (not checkDbl(scale,    "scale"))     return false;
    if (not checkStr(params,   "params"))    return false;
    if (not checkStr(model,    "model"))     return false;
    if (not checkStr(biorth,   "biorth"))    return false;

    // Arrays
    //
    file.getDataSet("xgrid").read(xgrid);

    // Allocate arrays
    //
    potl_array.resize(mmax+1);
    dens_array.resize(mmax+1);
    dpot_array.resize(mmax+1);
    rot_matrix.resize(mmax+1);

    if (cmap) {
      xmin = (rmin/scale - 1.0)/(rmin/scale + 1.0);
      xmax = (rmax/scale - 1.0)/(rmax/scale + 1.0);
      dxi = (xmax-xmin)/(numr-1);
    } else {
      xmin = rmin;
      xmax = rmax;
      dxi = (xmax-xmin)/(numr-1);
    }

    // Read from H5
    //
    for (int m=0; m<=mmax; m++) {
      std::ostringstream sout;
      sout << m;
      auto harmonic = file.getGroup(sout.str());
      
      harmonic.getDataSet("potl").read(potl_array[m]);
      harmonic.getDataSet("dens").read(dens_array[m]);
      harmonic.getDataSet("dpot").read(dpot_array[m]);
      harmonic.getDataSet("rot" ).read(rot_matrix[m]);
    }

  } catch (HighFive::Exception& err) {
    if (myid==0) std::cerr << "---- " << err.what() << std::endl;
    return false;
  }
    
  if (myid==0) std::cout << "---- EmpCyl2d::ReadH5Cache: "
			 << "read 2d basis cache <" << cache_name_2d << ">"
			 << std::endl;

  return true;
}


std::tuple<int, int, double, double> EmpCyl2d::linear_interp(double r)
{
  auto it = std::lower_bound(xgrid.begin(), xgrid.end(), r);

  int lo, hi;

  if (it == xgrid.begin()) {
    lo = 0;
    hi = 1;
  } else if (it == xgrid.end()) {
    hi = xgrid.size()-1;
    lo = hi - 1;
  } else {
    hi = std::distance(xgrid.begin(), it);
    lo = hi - 1;
  }

  double A = (xgrid[hi] - r)/(xgrid[hi] - xgrid[lo]);
  double B = (r - xgrid[lo])/(xgrid[hi] - xgrid[lo]);

  return {lo, hi, A, B};
}

void EmpCyl2d::checkMN(int& M, int& N, const std::string& member)
{
  if (M<0) M = 0;
  if (N<0) N = 0;
  if (M>mmax)  throw std::runtime_error("EmpCyl2d::"+member+": M too large");
  if (N>=nmax) throw std::runtime_error("EmpCyl2d::"+member+": N too large");
}

double EmpCyl2d::get_potl(double r, int M, int N)
{
  checkMN(M, N, "get_potl");

  if (basis_test) {
    return basis->potl(M, N, r)/sqrt(basis->norm(N, M));
  }

  int lo, hi;
  double A, B;
  std::tie(lo, hi, A, B) = linear_interp(r);

  return A*potl_array[M](lo, N) + B*potl_array[M](hi, N);
}

double EmpCyl2d::get_dens(double r, int M, int N)
{
  checkMN(M, N, "get_dens");

  if (basis_test) {
    return basis->dens(M, N, r)/sqrt(basis->norm(N, M));
  }

  int lo, hi;
  double A, B;
  std::tie(lo, hi, A, B) = linear_interp(r);

  return A*dens_array[M](lo, N) + B*dens_array[M](hi, N);
}

double EmpCyl2d::get_dpot(double r, int M, int N)
{
  checkMN(M, N, "get_dpot");

  if (basis_test) {
    return basis->dpot(M, N, r)/sqrt(basis->norm(N, M));
  }

  int lo, hi;
  double A, B;
  std::tie(lo, hi, A, B) = linear_interp(r);

  return A*dpot_array[M](lo, N) + B*dpot_array[M](hi, N);
}

void EmpCyl2d::get_pot(Eigen::MatrixXd& mat, double r)
{
  int lo, hi;			// Get the linear interp
  double A, B;
  std::tie(lo, hi, A, B) = linear_interp(r);

  mat.resize(mmax+1, nmax);	// Resize the return array

  for (int m=0; m<=mmax; m++) {	// Pack the array
    for (int j=0; j<nmax; j++)
      mat(m, j) = A*potl_array[m](lo, j) + B*potl_array[m](hi, j);
  }
}

void EmpCyl2d::get_dens(Eigen::MatrixXd& mat, double r)
{
  int lo, hi;
  double A, B;
  std::tie(lo, hi, A, B) = linear_interp(r);

  mat.resize(mmax+1, nmax);

  for (int m=0; m<=mmax; m++) {
    for (int j=0; j<nmax; j++)
      mat(m, j) = A*dens_array[m](lo, j) + B*dens_array[m](hi, j);
  }
}

void EmpCyl2d::get_force(Eigen::MatrixXd& mat, double r)
{
  int lo, hi;
  double A, B;
  std::tie(lo, hi, A, B) = linear_interp(r);

  mat.resize(mmax+1, nmax);

  for (int m=0; m<=mmax; m++) {
    for (int j=0; j<nmax; j++)
      mat(m, j) = A*dpot_array[m](lo, j) + B*dpot_array[m](hi, j);
  }
}

void EmpCyl2d::checkCoefs()
{
  if (myid) return;

  Mapping  map(scale, cmap);
  auto     disk = createModel();
  LegeQuad lw(knots);

  Eigen::VectorXd coefs(nmax), coef0(nmax);
  coefs.setZero();
  coef0.setZero();
  
  double ximin = map.r_to_xi(rmin);
  double ximax = map.r_to_xi(rmax);

  for (int k=0; k<knots; k++) {
    double xx  = ximin + (ximax - ximin)*lw.knot(k);
    double rr  = map.xi_to_r(xx);
    double fac = lw.weight(k) * rr / map.d_xi_to_r(xx) * (ximax - ximin);
    
      for (int j=0; j<nmax; j++) {
	coefs(j) += fac * disk->dens(rr) * get_potl(rr, 0, j);
	coef0(j) += fac * disk->dens(rr) * basis->potl(0, j, rr) / sqrt(basis->norm(j, 0));
      }
  }

  std::cout << std::endl << "Coefficients" << std::endl;

  for (int j=0; j<nmax; j++) {
    std::cout << std::setw( 6) << j
	      << std::setw(16) << coefs(j)
	      << std::setw(16) << coef0(j)
	      << std::endl;
  }

  std::cout << std::endl << "Reconstruction" << std::endl;

  double dx = (ximax - ximin)/(numr - 1);

  for (int k=0; k<numr; k++) {
    double xx  = ximin + dx*k;
    double rr  = map.xi_to_r(xx);
    double yy = 0.0, zz = 0.0, uu = 0.0, vv = 0.0;
    for (int j=0; j<nmax; j++) {
      uu += coefs(j) * get_potl(rr, 0, j);
      vv += coefs(j) * get_dens(rr, 0, j);
      yy += coef0(j) * basis->potl(0, j, rr) / sqrt(basis->norm(j, 0));
      zz += coef0(j) * basis->dens(0, j, rr) / sqrt(basis->norm(j, 0));
    }

    std::cout << std::setw(16) << rr
	      << std::setw(16) << uu
	      << std::setw(16) << vv
	      << std::setw(16) << yy
	      << std::setw(16) << zz
	      << std::endl;
  }
}
