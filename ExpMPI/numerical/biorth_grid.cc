// This may look like C code, but it is really -*- C++ -*-

//  Routines for computing biorthonormal pairs based on grid
//
//  MDW 11/13/91 [based on "findwake" by MDW: 12/26/1987]



static const char rcsid[] = "$Id$";

#include <string>
#include <iostream.h>
#include <math.h>
#include <Vector.h>
#include <interp.h>
#include "biorth.h"

double factrl(int n);
double plgndr(int l, int m, double x);

BiorthGrid::BiorthGrid(AxiSymBiorth& T, double RMIN, double RMAX, 
		       int NMAX, int LMAX, int RNUM)
{
  BiorthID = "BiorthGrid(" + T.BiorthID + ")";
  dof = T.get_dof();
  
  t = &T;

  rmin = RMIN;
  rmax = RMAX;
  lmax = LMAX;
  nmax = NMAX;
  rnum = RNUM;

  xmin = T.r_to_rb(rmin);
  xmax = T.r_to_rb(rmax);

  potl_grid  = new Matrix[lmax+1];
  potl_grid2 = new Matrix[lmax+1];
  dens_grid  = new Matrix[lmax+1];
  dens_grid2 = new Matrix[lmax+1];

  x_grid.setsize(1, rnum);

  krnl_grid.setsize(1, nmax, 0, lmax);
  norm_grid.setsize(0, nmax-1, 0, lmax);

  double del = (xmax - xmin)/(double)(rnum-1);
  for (int ir=0; ir<rnum; ir++) x_grid[ir+1] = xmin + del*ir;

  for (int l=0; l<=lmax; l++) {

    potl_grid[l].setsize(1, nmax, 1, rnum);
    potl_grid2[l].setsize(1, nmax, 1, rnum);
    dens_grid[l].setsize(1, nmax, 1, rnum);
    dens_grid2[l].setsize(1, nmax, 1, rnum);

    for (int n=1; n<=nmax; n++) {

      krnl_grid[n][l] = T.krnl(n, l);
//      norm_grid[n-1][l] = T.norm(n-1, l);
//      MDW 1/18/95: get_dens and get_potl are normalized by 1/sqrt(T.norm)
      norm_grid[n-1][l] = sqrt(T.norm(n-1, l));

      for (int ir=1; ir<=rnum; ir++) {
	potl_grid[l][n][ir] = T.potl(n, l, x_grid[ir]);
	dens_grid[l][n][ir] = T.dens(n, l, x_grid[ir]);
      }

      Spline(x_grid, potl_grid[l][n], -1.0e30, -1.0e30, potl_grid2[l][n]);
      Spline(x_grid, dens_grid[l][n], -1.0e30, -1.0e30, dens_grid2[l][n]);
    }
  }

}


double BiorthGrid::potl(int n, int l, double x)
{
  double ans;

  if (x > xmax || x < xmin) return 0.0;

  Splint1(x_grid, potl_grid[l][n], potl_grid2[l][n], x, ans, 1);

  return ans;
}

double BiorthGrid::dens(int n, int l, double x)
{
  double ans;

  if (x > xmax || x < xmin) return 0.0;

  Splint1(x_grid, dens_grid[l][n], dens_grid2[l][n], x, ans, 1);

  return ans;
}



double AxiSymBiorth::get_dens(double r, int l, const Vector& coef)
{
  double accum=0.0;

  for (int n=coef.getlow(); n<=coef.gethigh(); n++)
    accum += coef[n]*dens(n, l, r_to_rb(r))/sqrt(norm(n-1, l));

  return accum;

}

double AxiSymBiorth::get_potl(double r, int l, const Vector& coef)
{
  double accum=0.0;

  for (int n=coef.getlow(); n<=coef.gethigh(); n++)
    accum += coef[n]*potl(n, l, r_to_rb(r))/sqrt(norm(n-1, l));

  return -accum;

}


double plgndr(int l, int m, double x)
{
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;
  
  if (m < 0 || m > l || fabs(x) > 1.0) {
    cerr << "Bad arguments in routine PLGNDR\n";
    exit(-1);
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

double factrl(int n)
{
  if (n<=1)
    return 1;
  else
    return n*factrl(n-1);
}

