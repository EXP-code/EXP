/*
  compute scale length according to a Plummer profile
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void get_potacc(double,int,double *,double *,double *);

#define NUM 100

void adjust_scale(int n)
{

  int i;
  double r=0.0,dr,p,dp,p0;

  if (n%nscale) return;

  if (myid == 0) {

    get_potacc(0.0,0,expcoef[0],&p0,&dp);
    dr = rmax/NUM;
    for (i=1; i<=NUM; r +=dr, i++) {
      get_potacc(r,0,expcoef[0],&p,&dp);
      if (p>0.5*p0) {
	scale = r;
	break;
      }
    }

    MPI_Bcast(&scale, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  }
  else {

    MPI_Bcast(&scale, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  }

}

