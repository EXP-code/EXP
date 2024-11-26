#include <limits>

#include <Eigen/Eigen>

#include "expand.H"

#include <gaussQ.H>
#include <CBDisk.H>
#include <interp.H>

const std::set<std::string>
CBDisk::valid_keys = {
  "mmax",
  "scale",
  "diskconf",
  "background"
};

CBDisk::CBDisk(Component* c0, const YAML::Node& conf, MixtureBasis* m) :
  PolarBasis(c0, conf, m)
{
  id = "CBDisk";
  is_flat = true;

  // Radial scale factor/softening radius
  //
  scale = 1.0;

  // Get initialization info
  //
  initialize();

  if (myid==0) {
    std::string sep("----    ");
    std::cout << "---- CBDisk parameters: "
	      << std::endl << sep << "scale="       << scale
	      << std::endl << sep << "lmax="        << Lmax
	      << std::endl << sep << "nmax="        << nmax
	      << std::endl << sep << "NO_L0="       << std::boolalpha << NO_M0
	      << std::endl << sep << "NO_L1="       << std::boolalpha << NO_M1
	      << std::endl << sep << "EVEN_M="      << std::boolalpha << EVEN_M
	      << std::endl << sep << "M0_ONLY="     << std::boolalpha << M0_only
	      << std::endl << sep << "selfgrav="    << std::boolalpha << self_consistent
	      << std::endl;
  }

  setup();

  // Potential, force, and density scaling
  //
  fac1 = pow(scale, -0.5);
  fac2 = pow(scale, -1.5);

  if (myid==0) orthoCheck();
}


void CBDisk::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["mmax"])      mmax       = conf["mmax"].as<int>();
    if (conf["Lmax"])      mmax       = conf["Lmax"].as<int>();
    if (conf["Mmax"])      mmax       = conf["Mmax"].as<int>();
    if (conf["scale"])     scale      = conf["scale"].as<double>();

    Lmax = Mmax = mmax;		// Override base-class values
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in CBDisk: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("CBDisk::initialize: error in parsing YAML");
  }

  // Set background model
  if (conf["background"]) setBackground();
}

CBDisk::~CBDisk(void)
{
  // NADA
}


void CBDisk::setBackground()
{
  try {

    YAML::Node Params = conf["background"];

    std::string name = Params["name"].as<std::string>();
    auto params      = Params["parameters"];
    
    // Convert ID string to lower case
    //
    std::transform(name.begin(), name.end(), name.begin(),
		   [](unsigned char c){ return std::tolower(c); });
    
    if (name.find("kuzmin") != std::string::npos) {
      if (myid==0) {
	std::cout << "---- CBDisk: Making a Kuzmin disk";
	if (params) std::cout << " with" << std::endl << Params["parameters"];
	std::cout << std::endl;
      }
      disk = std::make_shared<EmpCyl2d::KuzminCyl>(params);
      return;
    }
      
    if (name.find("mestel") != std::string::npos) {
      if (myid==0) {
	std::cout << "---- CBDisk: Making a finite Mestel disk";
	if (params) std::cout << " with" << std::endl << Params["parameters"];
	std::cout << std::endl;
      }
      disk = std::make_shared<EmpCyl2d::MestelCyl>(params);
      return;
    }
      
    if (name.find("zang") != std::string::npos) {
      if (myid==0) {
	std::cout << "---- CBDisk: Making a double-tapered Zang";
	YAML::Emitter out;
	out << YAML::Flow << params;
	if (params) std::cout << " with " << out.c_str();
	std::cout << std::endl;
      }
      disk = std::make_shared<EmpCyl2d::ZangCyl>(params);
      return;
    }
      
    if (name.find("expon") != std::string::npos) {
      if (myid==0) {
	std::cout << "---- CBDisk: Making an Exponential disk";
	if (params) std::cout << " with" << std::endl << params;
	std::cout << std::endl;
      }
      disk = std::make_shared<EmpCyl2d::ExponCyl>(params);
      return;
    }
      
    // Default if nothing else matches
    if (myid==0)
      std::cout << "---- CBDisk: Making an Exponential disk [Default]"
		<< std::endl;
    disk = std::make_shared<EmpCyl2d::ExponCyl>(params);
    return;
  }
  catch (YAML::Exception & error) {
    if (myid==0)
      std::cout << "Error parsing parameters in CBDisk::setBackground: "
		<< error.what() << std::endl
		<< std::string(60, '-') << std::endl
		<< "Config node"        << std::endl
		<< std::string(60, '-') << std::endl
		<< conf["background"]   << std::endl
		<< std::string(60, '-') << std::endl;
    throw std::runtime_error("CBDisk::setBackground: error parsing YAML config");
  }
}


void CBDisk::get_dpotl(double r, double z,
		       Eigen::MatrixXd& pot,
		       Eigen::MatrixXd& dpr,
		       Eigen::MatrixXd& dpz, int tid)
{
  pot.resize(mmax+1, nmax);	// Resize the return arrays
  dpr.resize(mmax+1, nmax);
  dpz.resize(mmax+1, nmax);

  dpz.setZero();
  
  double R = r/scale;

  for (int m=0; m<=mmax; m++) {	// Pack the array
    for (int j=0; j<nmax; j++) {
      pot(m, j) = potl(j, m, R) * fac1;
      dpr(m, j) = -dpot(j, m, R) * fac1 / scale;
    }
  }
}

void CBDisk::get_potl(double r, double z, Eigen::MatrixXd& p, int tid)
{
  p.resize(mmax+1, nmax);	// Resize the return array

  double R = r/scale;

  for (int m=0; m<=mmax; m++) {	// Pack the array
    for (int j=0; j<nmax; j++)
      p(m, j) = potl(j, m, R) * fac1;
  }
}

void CBDisk::get_dens(double r, double z, Eigen::MatrixXd& p, int tid)
{
  p.resize(mmax+1, nmax);	// Resize the return array

  double R = r/scale;

  for (int m=0; m<=mmax; m++) {	// Pack the array
    for (int j=0; j<nmax; j++)
      p(m, j) = dens(j, m, R) * fac2;
  }
}

void CBDisk::get_potl_dens(double r, double z,
			   Eigen::MatrixXd& p, Eigen::MatrixXd& d,
			   int tid)
{
  p.resize(mmax+1, nmax);	// Resize the return arrays
  d.resize(mmax+1, nmax);

  double R = r/scale;

  for (int m=0; m<=mmax; m++) {	// Pack the array
    for (int j=0; j<nmax; j++) {
      p(m, j) = potl(j, m, R) * fac1;
      d(m, j) = dens(j, m, R) * fac2;
    }
  }
}

//  Routines for computing biorthonormal pairs based on 
//  Clutton-Brock's 2-dimensional series
//
double CBDisk::phif(const int n, const int m, const double r)
{
  // By recurrance relation
  //
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double cur = sqrt(fac);

  for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);

  if (n==0) return cur;
  
  double curl1 = cur;
  double curl2;
  
  fac *= r2 - 1.0;
  cur *= fac*(2*m+1);

  if (n==1) return cur;

  for (int nn=2; nn<=n; nn++) {
    curl2 = curl1;
    curl1 = cur;
    cur   = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
      (1.0 + (double)(2*m-1)/nn)*curl2;
  }

  return cur;
}

double CBDisk::potl(const int n, const int m, const double r)
{
  return pow(r, m) * phif(n, m, r)/sqrt(norm(n, m));
}

double CBDisk::dpot(const int n, const int m, const double r)
{
  double ret = dphi(n, m, r);
  if (m) ret = (phif(n, m, r)*m/r + ret) * pow(r, m);
  return ret/sqrt(norm(n, m));
}

double CBDisk::dphi(const int n, const int m, const double r)
{
  double ret = phif(n, m+1, r);
  if (n>0) ret -= 2.0*phif(n-1, m+1, r);
  if (n>1) ret += phif(n-2, m+1, r);
  return -r*ret;
}


// By recurrance relation
//
void CBDisk::potl(const int m, const double r, Eigen::VectorXd& a)
{
  a.resize(nmax);

  double pfac = pow(r, m);
  
  double r2  = r*r;
  double fac = 1.0/(1.0 + r2);
  double cur = sqrt(fac);

  for (int mm=1; mm<=m; mm++) cur *= fac*(2*mm - 1);

  a(0) = pfac*cur;

  if (nmax>0) {

    double curl1 = cur;
    double curl2;
  
    fac *= r2 - 1.0;
    cur *= fac*(2*m+1);

    a(1) = pfac*cur;

    for (int nn=2; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
	(1.0 + (double)(2*m-1)/nn)*curl2;
      a(nn) = pfac*cur;
    }
  }
  
  for (int n=0; n<nmax; n++) a(n) /= sqrt(norm(n, m));

  return;
}



double CBDisk::dens(int n, int m, double r)
{
  double f = 0.5/sqrt(norm(n, m))/M_PI;

  if (n>=2) 
    return f*pow(r, m)*(phif(n, m+1, r) - phif(n-2, m+1, r));
  else
    return f*pow(r, m)*phif(n, m+1, r);
}


void CBDisk::dens(int mm, double r, Eigen::VectorXd& a)
{
  a.resize(nmax);
  
  double pfac = pow(r, (double)mm+1.0e-20);
  
  int m = mm + 1;
  double r2 = r*r;
  double fac = 1.0/(1.0 + r2);
  double cur = sqrt(fac);

  for (int M=1; M<=m; M++) cur *= fac*(2*M - 1);

  a(0) = pfac*cur;

  if (nmax>0) {
  
    double curl1 = cur;
    double curl2;
  
    fac *= r2 - 1.0;
    cur *= fac*(2*m+1);

    a(1) = pfac*cur;

    for (int nn=2; nn<nmax; nn++) {
      curl2 = curl1;
      curl1 = cur;
      cur = (2.0 + (double)(2*m-1)/nn)*fac*curl1 - 
	(1.0 + (double)(2*m-1)/nn)*curl2;
      a(nn) = pfac*cur;
    }
    
    for (int nn=nmax-1; nn>1; nn--)
      a(nn) -= a(nn-2);
  }
  
  for (int n=0; n<nmax; n++) a(n) *= 0.5/sqrt(norm(n, mm))/M_PI;

  return;
}


double CBDisk::norm(int n, int m)
{
  double ans = 1.0;
  
  for (int i=n+1; i<=n+2*m; i++) ans *= i;

  return pow(0.5, 2*m+1)*ans;
}

void CBDisk::orthoCheck()
{
  const std::string hh("------    ");
  const double tol = 1.0e-4;
  const int num = 4000;
  LegeQuad lq(num);

  std::vector<Eigen::MatrixXd> ortho(mmax+1);
  Eigen::VectorXd worst(mmax+1);
  Eigen::MatrixXd vpot(mmax+1, nmax), vden(mmax+1, nmax);

  double Rmax = scale*100.0;

  worst.setZero();

  // Allocate result matrices
  //
  for (auto & v : ortho) {
    v.resize(nmax, nmax);
    v.setZero();
  }

  for (int i=0; i<num; i++) {
    double r = lq.knot(i) * Rmax, fac = lq.weight(i) * Rmax;

    get_potl(r, 0.0, vpot, 0);
    get_dens(r, 0.0, vden, 0);

    for (int m=0; m<=mmax; m++) {
      for (int j=0; j<nmax; j++) {
	for (int k=0; k<nmax; k++) {
	  ortho[m](j, k) += fac * vpot(m, j) * vden(m, k) * r * 2.0*M_PI;
	}
      }
    }
  }

  for (int m=0; m<=mmax; m++) {
    for (int j=0; j<nmax; j++) {
      for (int k=0; k<nmax; k++) {
	double test = ortho[m](j, k);
	if (j==k) test -= 1.0;
	if (fabs(test) > worst(m)) worst(m) = fabs(test);
      }
    }

    if (worst(m) > tol) {
      std::ostringstream sout;
      sout << "CBDisk.ortho." << m ;
      std::ofstream out(sout.str());
      if (out) out << ortho[m] << std::endl;
    }

  }

  std::cout << hh << std::endl
	    << hh << "Orthogonality check" << std::endl
	    << hh << std::endl
	    << hh << std::setw(6) << "m" << std::setw(16) << "worst" << std::endl
	    << hh << std::setw(6) << "-" << std::setw(16) << "-----" << std::endl;

  for (int m=0; m<=mmax; m++)
    std::cout << hh << std::setw(6) << m << std::setw(16) << worst(m) << std::endl;
  std::cout << std::endl;

  const int numr = 10;
  double dR = 3.0*scale/numr, dH = 0.01*scale;

  std::ofstream out("CBDisk.ortho.gradient");
  for (int m=0; m<=mmax; m++) {
    for (int i=0; i<numr; i++) {
      double r = dR*(0.5 + i);
      out << std::setw(4) << m << std::setw(16) << r;

      for (int j=0; j<std::min(nmax, 4); j++) {
	double pp = potl(j, m, (r+dH)/scale), pm = potl(j, m, (r-dH)/scale);
	out << std::setw(16) << (pp - pm)*0.5/dH/scale
	    << std::setw(16) << dpot(j, m, r/scale)/scale/scale;
      }
      out << std::endl;
    }
  }
}
