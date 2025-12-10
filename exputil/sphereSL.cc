#include <stdlib.h>
#include <iomanip>
#include "sphereSL.H"


int  SphereSL::mpi   = 0;
bool SphereSL::cache = true;

/**
   This is now a wrapper class for SLGridSph
   (why duplicate all of the effort and give myself a maintenance headache?)
*/
SphereSL::SphereSL(int LMAX, int NMAX, int NUMR,
		   double RMIN, double RMAX, double SCALE, 
		   std::shared_ptr<SphericalModelTable> mod)
{
  BiorthID = "SphereSL";
  dof = 3;

  lmax = LMAX;
  nmax = NMAX;
  numr = NUMR;

  rmin = RMIN;
  rmax = RMAX;
  scale = SCALE;

  SLGridSph::mpi = mpi;
  int cmap = 0;

  slgrid = new SLGridSph(mod, lmax, nmax, numr, rmin, rmax, cache, cmap, scale);

}

SphereSL::~SphereSL(void)
{
  delete slgrid;
}

				// Members

double SphereSL::r_to_rb(const double r)
{
  return slgrid->r_to_xi(r);
}
    
double SphereSL::rb_to_r(const double rb)
{
  return slgrid->xi_to_r(rb);
}

double SphereSL::d_r_to_rb(const double r)
{
  return 1.0/slgrid->d_xi_to_r(slgrid->r_to_xi(r));
}


double SphereSL::d_rb_to_r(const double rb)
{
  return slgrid->d_xi_to_r(rb);
}

double SphereSL::potl(const int n, const int l, const double xx)
{
  return slgrid->get_pot(xx, l, n);
}


double SphereSL::dens(const int n, const int l, const double xx)
{
  return slgrid->get_dens(xx, l, n);
}

double SphereSL::get_potl(const double r, const int l,
			  const Eigen::VectorXd& coef)
{
  double x = r_to_rb(r);
  double ans=0.0;
  Eigen::VectorXd val(nmax);

  slgrid->get_pot(val, x, l);
  for (int n=1; n<=nmax; n++) ans += val[n]*coef[n];

  return ans;
}


double SphereSL::get_dens(const double r, const int l, const Eigen::VectorXd& coef)
{
  double x = r_to_rb(r);

  double ans=0.0;
  Eigen::VectorXd val(1, nmax);

  slgrid->get_dens(val, x, l);
  for (int n=1; n<=nmax; n++) ans += val[n]*coef[n];

  return ans;
}

void SphereSL::potl(const int nn, const int l, const double x,
		    Eigen::VectorXd& t)
{
  t.resize(nmax);
  slgrid->get_pot(t, x, l);
}


void SphereSL::dens(const int nn, const int l, const double x,
		    Eigen::VectorXd& t)
{
  t.resize(nmax);
  slgrid->get_dens(t, x, l);
}


