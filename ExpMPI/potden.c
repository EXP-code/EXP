
#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

				/* Return radial wave function for potential */
double potli(double rr, int l, int n)
{
  double r,temp,d1,d2,r_to_rh(double r),r_to_rq(double r);

  if (c_brock)
    r = r_to_rh(rr);
  else if (hernq)
    r = r_to_rq(rr);
  else
    r = rr;

  splint(r_grid,potl_grid[l].rw[n],potl_grid[l].rw2[n],RNUM,r,&temp,&d1,&d2);

  return temp;
}

				/* Return pointers to radial wave function
				   and its derivative for potential */
void potlgri(double rr, int l, int n, double *p, double *dp)
{
  double r,d2;
  double r_to_rh(double r),d_r_to_rh(double r);
  double r_to_rq(double r),d_r_to_rq(double r);

  if (c_brock)
    r = r_to_rh(rr);
  else if (hernq)
    r = r_to_rq(rr);
  else
    r = rr;

  splint(r_grid,potl_grid[l].rw[n],potl_grid[l].rw2[n],RNUM,r,p,dp,&d2);

  if (c_brock)
    *dp *= d_r_to_rh(rr);
  if (hernq)
    *dp *= d_r_to_rq(rr);
}


				/* Return radial wave function for density   */
				/* [Don't forget the 4*PI dimensional factor */
double densi(double rr, int l, int n)
{
  double r,temp,d1,d2,r_to_rh(double r),r_to_rq(double r);

  if (c_brock)
    r = r_to_rh(rr);
  else if (hernq)
    r = r_to_rq(rr);
  else
    r = rr;

  splint(r_grid,dens_grid[l].rw[n],dens_grid[l].rw2[n],RNUM,r,&temp,&d1,&d2);

  return temp;
}



				/* Return pointers for potential and 
				   derivatives for given coef vector */

void get_potacc(double r, int l, double *coef, double *p, double *dp)
{
  int n;
  double accum1=0.0,accum2=0.0,p1,dp1;

  for (n=1; n<=potl_grid[l].nmax; n++) {
    potlgri(r,l,n,&p1,&dp1);
    accum1 += coef[n]*p1;
    accum2 += coef[n]*dp1;
  }

  *p = -accum1;
  *dp = -accum2;
}

