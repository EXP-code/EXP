
/*
  Compute acceleration and potential for on particles due
  to point masses
*/

static char rcsid[] = "$Id$";

#include "expand.h"

void make_pointmass_list()
{
  MPI_Status status;
  int i, *p;
  static int firstime = 1;


  pmnum = 0;			/* Count point masses */
  for (i=1; i<=nbodies; i++)
    if (component[i] == 0) pmnum++;

#ifdef DEBUG
  fprintf(stderr, "Process %d: nbodies=%d  pmnum=%d\n", myid, nbodies, pmnum);
#endif
				/* Delete point mass list */
  if (!firstime) free(pmlist);

				/* If there are point masses, make list */
  if (pmnum>0) {
    pmlist = (int *)malloc(pmnum*sizeof(int));
    firstime = 0;
    p = pmlist;
    for (i=1; i<=nbodies; i++)
      if (component[i] == 0) *(p++) = i;    
  }


				/* Broadcast total number of point masses */
  MPI_Allreduce(&pmnum, &pmnum0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  /*
  MPI_Reduce(&pmnum, &pmnum0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==0) {
    for (i=1; i<numprocs; i++)
      MPI_Send(&pmnum0, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
  } else {
      MPI_Recv(&pmnum0, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  }
  */

}

void get_acceleration_and_potential_pointmass()
{
  double fac, ffac, xx, yy, zz, mm;
  int i, j, k, l;
  int pmnum1;

  for (l=0; l<numprocs; l++) {
    
    if (l == myid) pmnum1 = pmnum;

				/* Send number of particles */
    MPI_Bcast(&pmnum1, 1, MPI_INT, l, MPI_COMM_WORLD);

				/* Skip if zero */
    for (j=0; j<pmnum1; j++) {

				/* Send mass and pos'ns */
      if (l == myid) {

	k = pmlist[j];

	mm = mass[k];
	xx = x[k];
	yy = y[k];
	zz = z[k];
      }

      MPI_Bcast(&mm, 1, MPI_DOUBLE, l, MPI_COMM_WORLD);
      MPI_Bcast(&xx, 1, MPI_DOUBLE, l, MPI_COMM_WORLD);
      MPI_Bcast(&yy, 1, MPI_DOUBLE, l, MPI_COMM_WORLD);
      MPI_Bcast(&zz, 1, MPI_DOUBLE, l, MPI_COMM_WORLD);


      for (i=1; i<=nbodies; i++) {
				/* No self force, please */
	if (l == myid && k == i) continue;

	if (freeze_particle(i)) continue;

	fac = pow(pmsoft*pmsoft + 
		  (x[i]-xx)*(x[i]-xx) +
		  (y[i]-yy)*(y[i]-yy) +
		  (z[i]-zz)*(z[i]-zz), -0.5);
		     
	ffac = -mm*fac*fac*fac;

	ax[i] += ffac * (x[i]-xx);
	ay[i] += ffac * (y[i]-yy);
	az[i] += ffac * (z[i]-zz);
	potext[i] += -mm*fac;
      }
    
    }

  }

}
