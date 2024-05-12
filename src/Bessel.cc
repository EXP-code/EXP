#include "expand.H"

#include <cmath>

#include <interp.H>
#include <Bessel.H>
#include <EXPmath.H>

const std::set<std::string>
Bessel::valid_keys = {
  "rnum"
};

void Bessel::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["rnum"])      RNUM       = conf["rnum"].as<int>();
    else                   RNUM       = 1000;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Sphere: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("Sphere::initialze: error parsing YAML");
  }
}

Bessel::Bessel(Component* c0, const YAML::Node& conf, MixtureBasis* m) : SphericalBasis(c0, conf, m)
{
  id = "BesselForce";
  initialize();

  // Initialize radial grids

  make_grid(0, rmax, Lmax, nmax);

  setup();
}

// Get potential functions by from table
void Bessel::get_dpotl(int lmax, int nmax, double r, 
		       Eigen::MatrixXd& p, Eigen::MatrixXd& dp, int tid)
{
  int klo = (int)( (r-r_grid[1])/r_grid_del ) + 1;
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;
  
  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;
  double aaa = -(3.0*a*a - 1.0)*r_grid_del/6.0;
  double bbb =  (3.0*b*b - 1.0)*r_grid_del/6.0;

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = a*potl_grid[l].rw(n, klo) + b*potl_grid[l].rw(n, khi) +
	aa*potl_grid[l].rw2(n, klo) + bb*potl_grid[l].rw2(n, khi);
      dp(l, n) = (-potl_grid[l].rw(n, klo)+potl_grid[l].rw(n, khi))/r_grid_del+
	aaa*potl_grid[l].rw2(n, klo) + bbb*potl_grid[l].rw2(n, khi);
    }
  }
}

void Bessel::get_potl(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  int klo = (int)( (r-r_grid[1])/r_grid_del ) + 1;
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;

  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = a*potl_grid[l].rw(n, klo) + b*potl_grid[l].rw(n, khi) +
	aa*potl_grid[l].rw2(n, klo) + bb*potl_grid[l].rw2(n, khi);
    }
  }
}

void Bessel::get_dens(int lmax, int nmax, double r, Eigen::MatrixXd& p, int tid)
{
  int klo = (int)( (r-r_grid[1])/r_grid_del ) + 1;
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;

  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = a*dens_grid[l].rw(n, klo) + b*dens_grid[l].rw(n, khi) +
	aa*dens_grid[l].rw2(n, klo) + bb*dens_grid[l].rw2(n, khi);
    }
  }
}


void Bessel::get_potl_dens(int lmax, int nmax, double r, 
			   Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid)
{
  int klo = (int)( (r-r_grid[1])/r_grid_del ) + 1;
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;

  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = a*potl_grid[l].rw(n, klo) + b*potl_grid[l].rw(n, khi) +
	aa*potl_grid[l].rw2(n, klo) + bb*potl_grid[l].rw2(n, khi);
      d(l, n) = a*dens_grid[l].rw(n, klo) + b*dens_grid[l].rw(n, khi) +
	aa*dens_grid[l].rw2(n, klo) + bb*dens_grid[l].rw2(n, khi);
    }
  }
}

double Bessel::dens(double r, int n)
{
  double alpha;

  if (n>p->n)
    throw GenericError("Routine dens() called with n out of bounds", __FILE__, __LINE__, 1001, true);

  alpha = p->a[n];
  return alpha*M_SQRT2/fabs(EXPmath::sph_bessel(p->l+1, alpha)) * pow(rmax,-2.5) *
    EXPmath::sph_bessel(p->l, alpha*r/rmax);
}

double Bessel::potl(double r, int n)
{
  double alpha;

  if (n>p->n)
    throw GenericError("Routine potl() called with n out of bounds", __FILE__, __LINE__, 1002, true);

  alpha = p->a[n];
  return M_SQRT2/fabs(alpha*EXPmath::sph_bessel(p->l+1,alpha)) * pow(rmax,-0.5) *
    EXPmath::sph_bessel(p->l,alpha*r/rmax);
}

void Bessel::make_grid(double rmin, double rmax, int lmax, int nmax)
{
  
  potl_grid.resize(lmax+1);
  dens_grid.resize(lmax+1);
  
  r_grid.resize(RNUM);

  r_grid_del = rmax/(double)(RNUM-1);
  double r = 0.0;
  for (int ir=0; ir<RNUM; ir++, r+=r_grid_del) r_grid[ir] = r;

  for (int l=0; l<=lmax; l++) {
    potl_grid[l].rw .resize(nmax, RNUM);
    potl_grid[l].rw2.resize(nmax, RNUM);
    dens_grid[l].rw .resize(nmax, RNUM);
    dens_grid[l].rw2.resize(nmax, RNUM);

    p = std::make_shared<Roots>(l, nmax);

    for (int n=0; n<nmax; n++) {
      r = 0.0;
      for (int ir=0; ir<RNUM; ir++, r+=r_grid_del) {
	potl_grid[l].rw(n, ir) = potl(r, n);
	dens_grid[l].rw(n, ir) = dens(r, n);
      }
      
      {
	Eigen::VectorXd X(potl_grid[l].rw.row(n));
	Eigen::VectorXd Y(X.size());
	Spline(r_grid, X, 1.0e30, 1.0e30, Y);
	potl_grid[l].rw2.row(n) = Y;
      }

      {
	Eigen::VectorXd X(dens_grid[l].rw.row(n));
	Eigen::VectorXd Y(X.size());
	Spline(r_grid, X, 1.0e30, 1.0e30, Y);
	dens_grid[l].rw2.row(n) = Y;
      }
    }

    potl_grid[l].nmax = nmax;
    dens_grid[l].nmax = nmax;
  }

  // check table

  for (int ir=0; ir<RNUM; ir++) assert(!std::isnan(r_grid[ir]));
  for (int l=0; l<=lmax; l++) {
    assert(!std::isnan(potl_grid[l].nmax));
    assert(!std::isnan(dens_grid[l].nmax));
    for (int n=0; n<nmax; n++) {
      for (int ir=0; ir<RNUM; ir++) {
	assert(!std::isnan(potl_grid[l].rw (n, ir)));
	assert(!std::isnan(potl_grid[l].rw2(n, ir)));
	assert(!std::isnan(dens_grid[l].rw (n, ir)));
	assert(!std::isnan(dens_grid[l].rw2(n, ir)));
      }
    }
  }
  
}

