/*****************************************************************************
 *  Description:
 *  -----------
 *
 *
 *  Routines for computing biorthonormal pairs based on j_n(x)
 *
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91 [based on "findwake" by MDW: 12/26/1987]
 *
 ***************************************************************************/

#include "expand.h"
#include <assert.h>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif


/*   Set up a structure for roots of a given order for biorthonormal 
     functions based on bessels */
struct ROOTS {
  int l;
  int n;
  double *a;
} *get_ln(int l, int n);

void rid_ln(struct ROOTS *p);
void set_radius(double R);
double dens(double r, int n, struct ROOTS *p);
double potl(double r, int n, struct ROOTS *p);

/*   Radius of system  */
static double RR;

void set_radius(double R)
{
  RR = R;
}

struct ROOTS *get_ln(int l, int n)
{
  struct ROOTS *p;
  double *a,*sbessjz();

  if ((p = (struct ROOTS *)malloc(sizeof(struct ROOTS))) == NULL)
    myerror("allocation failure in get_ln()");

  a = sbessjz(l-1,n);

  p->a = a;
  p->l = l;
  p->n = n;

  return p;
}

void rid_ln(struct ROOTS *p)
{
  free_dvector(p->a, 1, p->n);
  free(p);
}


/* given the appropriate ROOTS structure, make dens and pot functions */
#define SQ2 1.4142135623730950488
double dens(double r, int n, struct ROOTS *p)
{
  double sbessj(),alpha;

  if (n>p->n)
    myerror("Routine dens() called with n out of bounds");

  alpha = p->a[n];
  return alpha*SQ2/fabs(sbessj(p->l,alpha)) * pow(RR,-2.5) *
    sbessj(p->l,alpha*r/RR);
}

double potl(double r, int n, struct ROOTS *p)
{
  double sbessj(),alpha;

  if (n>p->n)
    myerror("Routine potl() called with n out of bounds");

  alpha = p->a[n];
  return SQ2/fabs(alpha*sbessj(p->l,alpha)) * pow(RR,-0.5) *
    sbessj(p->l,alpha*r/RR);
}

void make_grid(double rmin, double rmax, int lmax, int nmax)
{
  double r;
  int l,n,ir;
  struct ROOTS *p;

  potl_grid = (struct RGRID *) malloc((unsigned) (lmax+1)*sizeof(struct RGRID));
  if (!potl_grid) {
    fprintf(stderr, "biorth_bes: problem allocating <potl_grid>\n");
    exit(-1);
  }
  dens_grid = (struct RGRID *) malloc((unsigned) (lmax+1)*sizeof(struct RGRID));
  if (!dens_grid) {
    fprintf(stderr, "biorth_bes: problem allocating <dens_grid>\n");
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
    p = get_ln(l, nmax);
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

  /* check table */

  for (ir=1; ir<=RNUM; ir++) assert(!isnan(r_grid[ir]));
  for (l=0; l<=lmax; l++) {
    assert(!isnan(potl_grid[l].nmax));
    assert(!isnan(dens_grid[l].nmax));
    for (n=1; n<=nmax; n++) {
      for (ir=1; ir<=RNUM; ir++) {
	assert(!isnan(potl_grid[l].rw[n][ir]));
	assert(!isnan(potl_grid[l].rw2[n][ir]));
	assert(!isnan(dens_grid[l].rw[n][ir]));
	assert(!isnan(dens_grid[l].rw2[n][ir]));
      }
    }
  }
  fprintf(stderr,"Process %d: biorth_bes table OK\n", myid);

}

