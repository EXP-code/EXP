#include "expand.h"

#include <Vector.h>

void parallel_gather_coefficients(Matrix& expcoef, Matrix& expcoef1,
				  Matrix*& cc, Matrix*& cc1,
				  int lmax)
{
  int Ldim, L0, loffset, moffset, l, m, n, nn;

#ifdef MPE_PROFILE
  MPE_Log_event(5, myid, "b_gather_c");
#endif

  if (dof==3) {
    Ldim = lmax*(lmax + 2) + 1;
    L0 = 0;
  }
  else {
    Ldim = 2*lmax + 1;
    L0 = lmax;
  }    

  if (myid == 0) {

    for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      for (m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  for (n=1; n<=nmax; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;

	    for (nn=n; nn<=nmax; nn++)
	      cc[loffset+moffset][n][nn] = 0.0;
	  }
	  moffset++;
	}
	else {
	  for (n=1; n<=nmax; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;
	    expcoef[loffset+moffset+1][n] = 0.0;

	    for (nn=n; nn<=nmax; nn++) {
	      cc[loffset+moffset][n][nn] = 0.0;
	      cc[loffset+moffset+1][n][nn] = 0.0;
	    }
	  }
	  moffset+=2;
	}
      }
    }
  }


  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    for (m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for (n=1; n<=nmax; n++)
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	moffset++;
      }
      else {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&expcoef1[loffset+moffset+1][1],
		   &expcoef[loffset+moffset+1][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	for (n=1; n<=nmax; n++) {
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&cc1[loffset+moffset+1][n][n],
		     &cc[loffset+moffset+1][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	moffset+=2;
      }
    }
  }

#ifdef MPE_PROFILE
  MPE_Log_event(6, myid, "e_gather_c");
#endif

}

void parallel_distribute_coefficients(Matrix& expcoef, int lmax)
{
  int Ldim, L0, loffset, moffset, l, m;

#ifdef MPE_PROFILE
  MPE_Log_event(7, myid, "b_distrib_c");
#endif

  if (dof==3) {
    Ldim = lmax*(lmax + 2) + 1;
    L0 = 0;
  }
  else {
    Ldim = 2*lmax + 1;
    L0 = lmax;
  }    

  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      for (m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  MPI_Bcast(&expcoef[loffset+moffset][1], nmax, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset++;
	}
	else {
	  MPI_Bcast(&expcoef[loffset+moffset][1], nmax, MPI_DOUBLE,
		     0, MPI_COMM_WORLD);
	  MPI_Bcast(&expcoef[loffset+moffset+1][1], nmax, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset+=2;
	}
      }
  }

#ifdef MPE_PROFILE
  MPE_Log_event(8, myid, "e_distrib_c");
#endif

}

