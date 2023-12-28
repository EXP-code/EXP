//  Routines for computing biorthonormal pairs based on grid
//
//  MDW 11/13/91 [based on "findwake" by MDW: 12/26/1987]


#include <string>
#include <iostream>
#include <cmath>

#include <interp.H>
#include <biorth.H>

double factrl(int n);
double plgndr(int l, int m, double x);

BiorthGrid::BiorthGrid(AxiSymBioPtr T, double RMIN, double RMAX, 
		       int NMAX, int LMAX, int RNUM)
{
  BiorthID = "BiorthGrid(" + T->BiorthID + ")";
  dof = T->get_dof();
  
  t = T;

  rmin = RMIN;
  rmax = RMAX;
  lmax = LMAX;
  nmax = NMAX;
  rnum = RNUM;

  xmin = T->r_to_rb(rmin);
  xmax = T->r_to_rb(rmax);

  potl_grid  .resize(lmax+1);
  potl_grid2 .resize(lmax+1);
  dens_grid  .resize(lmax+1);
  dens_grid2 .resize(lmax+1);

  x_grid.resize(rnum);

  krnl_grid.resize(nmax, lmax);
  norm_grid.resize(nmax, lmax);

  double del = (xmax - xmin)/(double)(rnum-1);
  for (int ir=0; ir<rnum; ir++) x_grid[ir] = xmin + del*ir;

  for (int l=0; l<=lmax; l++) {

    potl_grid[l] .resize(nmax, rnum);
    potl_grid2[l].resize(nmax, rnum);
    dens_grid[l] .resize(nmax, rnum);
    dens_grid2[l].resize(nmax, rnum);

    for (int n=0; n<nmax; n++) {

      krnl_grid(n, l) = T->krnl(n, l);
      norm_grid(n, l) = sqrt(T->norm(n, l));

      for (int ir=0; ir<rnum; ir++) {
	potl_grid[l](n, ir) = T->potl(n, l, x_grid[ir]);
	dens_grid[l](n, ir) = T->dens(n, l, x_grid[ir]);
      }

      Eigen::VectorXd work(potl_grid2[l].row(n).size());

      Spline(x_grid, potl_grid[l].row(n), -1.0e30, -1.0e30, work);
      potl_grid2[l].row(n) = work;

      Spline(x_grid, dens_grid[l].row(n), -1.0e30, -1.0e30, work);
      dens_grid2[l].row(n) = work;
    }
  }

}


double BiorthGrid::potl(int n, int l, double x)
{
  double ans;

  if (x > xmax || x < xmin) return 0.0;

  Splint1(x_grid, potl_grid[l].row(n), potl_grid2[l].row(n), x, ans, 1);

  return ans;
}

double BiorthGrid::dens(int n, int l, double x)
{
  double ans;

  if (x > xmax || x < xmin) return 0.0;

  Splint1(x_grid, dens_grid[l].row(n), dens_grid2[l].row(n), x, ans, 1);

  return ans;
}



double AxiSymBiorth::get_dens(double r, int l, const Eigen::VectorXd& coef)
{
  double accum=0.0;

  for (int n=0; n<coef.size(); n++)
    accum += coef[n]*dens(n, l, r_to_rb(r))/sqrt(norm(n, l));

  return accum;

}

double AxiSymBiorth::get_potl(double r, int l, const Eigen::VectorXd& coef)
{
  double accum=0.0;

  for (int n=0; n<coef.size(); n++)
    accum += coef[n]*potl(n, l, r_to_rb(r))/sqrt(norm(n, l));

  return -accum;

}


double plgndr(int l, int m, double x)
{
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;
  
  if (m < 0 || m > l || fabs(x) > 1.0) {
    std::ostringstream sout;
    sout << "plgndr: bad arguments in routine l=" << l
	 << " m=" << m << " x=" << x;
    throw std::runtime_error(sout.str());
  }

  pmm=1.0;

  if (m > 0) {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=(m+2);ll<=l;ll++) {
	pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	pmm=pmmp1;
	pmmp1=pll;
      }
      return pll;
    }
  }
}


double dplgndr(int l, int m, double x)
{
  if (m < 0 || m > l || fabs(x) > 1.0) {
    std::ostringstream sout;
    sout << "dplgndr: bad arguments in routine l=" << l
	 << " m=" << m << " x=" << x;
    throw std::runtime_error(sout.str());
  }

  if (l==0 && m==0) return 0.0;

  double somx2 = 1.0/(x*x - 1.0);

  if (m<l)
    return somx2*(x*l*plgndr(l, m, x) - plgndr(l-1, m, x)*(l+m));
  else
    return somx2*x*l*plgndr(l, l, x);
}


double factrl(int n)
{
  if (n<=1)
    return 1;
  else
    return n*factrl(n-1);
}

