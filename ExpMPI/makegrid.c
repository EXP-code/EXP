#include "expand.h"
#include <assert.h>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

/******Definitions for routines <biorth_bes.c>******/

/*   Set up a structure for roots of a given order for biorthonormal 
     functions based on bessels */
struct ROOTS {
  int l;
  int n;
  double *a;
} *get_ln(int l, int n);


void rid_ln(struct ROOTS *p);
void set_radius(double R);
double dens(double r, int n, struct ROOTS *p),potl(double r, int n, struct ROOTS *p);

void make_grid(double rmin, double rmax, int lmax, int nmax)
{
  double r;
  int l,n,ir;
  struct ROOTS *p;

  potl_grid = (struct RGRID *) malloc((unsigned) (lmax+1)*sizeof(struct RGRID));
  if (!potl_grid) {
    fprintf(stderr, "bessacp: problem allocating <potl_grid>\n");
    exit(-1);
  }
  dens_grid = (struct RGRID *) malloc((unsigned) (lmax+1)*sizeof(struct RGRID));
  if (!dens_grid) {
    fprintf(stderr, "bessacp: problem allocating <dens_grid>\n");
    exit(-1);
  }

  r_grid = dvector(1,RNUM);

  r_grid_del = rmax/(double)(RNUM-1);
  for (ir=1, r=0.0; ir<=RNUM; ir++, r+=r_grid_del)
    r_grid[ir] = r;

  for (l=0; l<=lmax; l++) {
    potl_grid[l].rw = dmatrix(1,nmax,1,RNUM);
    potl_grid[l].rw2 = dmatrix(1,nmax,1,RNUM);
    dens_grid[l].rw = dmatrix(1,nmax,1,RNUM);
    dens_grid[l].rw2 = dmatrix(1,nmax,1,RNUM);
    p = get_ln(l,nmax);
    for (n=1; n<=nmax; n++) {
      for (ir=1, r=0.0; ir<=RNUM; ir++, r+=r_grid_del) {
	potl_grid[l].rw[n][ir] = potl(r,n,p);
	dens_grid[l].rw[n][ir] = dens(r,n,p);
      }
      
      spline(r_grid,potl_grid[l].rw[n],RNUM,1.0e30,1.0e30,
	     potl_grid[l].rw2[n]);
      spline(r_grid,dens_grid[l].rw[n],RNUM,1.0e30,1.0e30,
	     dens_grid[l].rw2[n]);
    }
    rid_ln(p);
    potl_grid[l].nmax = nmax;
    dens_grid[l].nmax = nmax;
  }
}

