#include <cassert>

#include <BiorthBess.H>
#include <EXPmath.H>
#include <EXPException.H>
#include <interp.H>
#include <gaussQ.H>
#include <cassert>

BiorthBess::BiorthBess(double rmax, int lmax, int nmax, int RNUM) :
  rmax(rmax), lmax(lmax), nmax(nmax), RNUM(RNUM)
{
  // Initialize radial grids
  make_grid();
}

// Get potential functions by from table
void BiorthBess::get_dpotl(double r, Eigen::MatrixXd& p)
{
  int klo = (int)( (r-r_grid[0])/r_grid_del );
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;
  
  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;
  double aaa = -(3.0*a*a - 1.0)*r_grid_del/6.0;
  double bbb =  (3.0*b*b - 1.0)*r_grid_del/6.0;

  p.resize(lmax+1, nmax);

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = (-potl_grid[l].rw(n, klo)+potl_grid[l].rw(n, khi))/r_grid_del+
	aaa*potl_grid[l].rw2(n, klo) + bbb*potl_grid[l].rw2(n, khi);
    }
  }
}

void BiorthBess::get_potl(double r, Eigen::MatrixXd& p)
{
  int klo = (int)( (r-r_grid[0])/r_grid_del );
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;

  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;

  p.resize(lmax+1, nmax);

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = a*potl_grid[l].rw(n, klo) + b*potl_grid[l].rw(n, khi) +
	aa*potl_grid[l].rw2(n, klo) + bb*potl_grid[l].rw2(n, khi);
    }
  }
}

void BiorthBess::get_dens(double r, Eigen::MatrixXd& p)
{
  int klo = (int)( (r-r_grid[0])/r_grid_del );
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;

  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;

  p.resize(lmax+1, nmax);

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = a*dens_grid[l].rw(n, klo) + b*dens_grid[l].rw(n, khi) +
	aa*dens_grid[l].rw2(n, klo) + bb*dens_grid[l].rw2(n, khi);
    }
  }
}


void BiorthBess::get_potl_dens(double r, Eigen::MatrixXd& p, Eigen::MatrixXd& d)
{
  int klo = (int)( (r-r_grid[0])/r_grid_del );
  if (klo < 0) klo = 0;
  if (klo > RNUM - 2) klo = RNUM - 2;
  int khi = klo + 1;

  double a = (r_grid[khi] - r)/r_grid_del;
  double b = (r - r_grid[klo])/r_grid_del;

  double aa = a*(a*a-1.0)*r_grid_del*r_grid_del/6.0;
  double bb = b*(b*b-1.0)*r_grid_del*r_grid_del/6.0;

  p.resize(lmax+1, nmax);
  d.resize(lmax+1, nmax);

  for (int l=0; l<=lmax; l++) {
    for (int n=0; n<nmax; n++) {
      p(l, n) = a*potl_grid[l].rw(n, klo) + b*potl_grid[l].rw(n, khi) +
	aa*potl_grid[l].rw2(n, klo) + bb*potl_grid[l].rw2(n, khi);
      d(l, n) = a*dens_grid[l].rw(n, klo) + b*dens_grid[l].rw(n, khi) +
	aa*dens_grid[l].rw2(n, klo) + bb*dens_grid[l].rw2(n, khi);
    }
  }
}

double BiorthBess::dens(double r, int n)
{
  double alpha;

  if (n>p->n)
    throw GenericError("Routine dens() called with n out of bounds", __FILE__, __LINE__, 1001, true);

  alpha = p->a[n];
  return alpha*M_SQRT2/fabs(EXPmath::sph_bessel(p->l+1, alpha)) * pow(rmax,-2.5) *
    EXPmath::sph_bessel(p->l, alpha*r/rmax);
}

double BiorthBess::potl(double r, int n)
{
  double alpha;

  if (n>p->n)
    throw GenericError("Routine potl() called with n out of bounds", __FILE__, __LINE__, 1002, true);

  alpha = p->a[n];
  return M_SQRT2/fabs(alpha*EXPmath::sph_bessel(p->l+1,alpha)) * pow(rmax,-0.5) *
    EXPmath::sph_bessel(p->l,alpha*r/rmax);
}

void BiorthBess::make_grid()
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

std::vector<Eigen::MatrixXd> BiorthBess::orthoCheck(int num)
{
  // Allocate storage
  //
  std::vector<Eigen::MatrixXd> one(lmax+1);
  for (auto & u : one) {
    u.resize(nmax, nmax);
    u.setZero();
  }

  // Number of knots
  //
  LegeQuad wk(num);
  
  // Radial range
  //
  double dr = rmax - rmin;

  Eigen::MatrixXd p(lmax+1, nmax), d(lmax+1, nmax);


  // Biorthogonal integral loop
  //
  for (int i=0; i<num; i++) {

    double r = rmin + dr*wk.knot(i);
    double w = dr*wk.weight(i) * r * r;

    // Evaluate basis at radius r
    //
    get_potl(r, p);
    get_dens(r, d);

    // Contribution to the integrand
    //
    for (int L=0, cnt=0; L<=lmax; L++) {
      for (int n1=0; n1<nmax; n1++) {
	for (int n2=0; n2<nmax; n2++, cnt++) {
	  one[L](n1, n2) += w * p(L, n1) * d(L, n2);
	}
      }
    }
  }

  return one;
}
