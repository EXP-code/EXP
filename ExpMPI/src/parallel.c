/*

  Parallel implementation routines

*/

/*
  Tags:

  10 = nbodies
  11 = restart
  12 = newrate
  18 = component
  19 = mass0
  20 = mass
  21 = x
  22 = y
  23 = z
  24 = vx
  25 = vy
  26 = vz
  27 = ax
  28 = ay
  29 = az
  30 = pot
  31 = potext
  32 = esave
  33 = mfp

*/

#include <stdio.h>
#include <math.h>
#include "expand.h"
#include "localmpi.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

				/* Round rate up unless storage */
				/* array will be overfilled     */
static double maxrate;

int round_up(double dnumb)
{
  int numb = dnumb + 1.0;
  if (numb >= nbodmax) numb = nbodmax;
  return numb;
}

void signal_failure(int fail)
{
  MPI_Bcast(&fail, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (fail) {
    MPI_Finalize();
    exit(-1);
  }
}


void setup_distribution(void)
{
  FILE *in, *out;
  int n;
  double norm=0.0;

  int fail = 0, fail0 = 0;

				/* Needed for both master and slaves */

  nbodies_index = (int *)malloc( numprocs * sizeof(int) );
  if (nbodies_index==NULL) {
    fprintf(stderr, "Process %d: Error inializing index\n", myid);
    fail = 1;
  }

  nbodies_table = (int *)malloc( numprocs * sizeof(int) );
  if (nbodies_table==NULL) {
    fprintf(stderr, "Process %d: Error inializing table\n", myid);
    fail = 1;
  }

  rates = (double *)malloc( numprocs * sizeof(double) );
  if (rates==NULL) {
    fprintf(stderr, "Process %d: Error inializing rates\n", myid);
    fail = 1;
  }

  if (myid == 0) {

    if ( (in = fopen("processor.rates", "r")) == NULL ) {
      fprintf(stderr, "setup: Error opening <processor.rates> . . . quit\n");
      fail = 1;
    }

    orates = (double *)malloc( numprocs * sizeof(double) );
    if (orates==NULL) {
      fprintf(stderr, "setup: Error inializing orates\n");
      fail = 1;
    }

    trates = (double *)malloc( numprocs * sizeof(double) );
    if (trates==NULL) {
      fprintf(stderr, "setup: Error inializing trates\n");
      fail = 1;
    }

    for (n=0; n<numprocs; n++) {
      if (fscanf(in, "%lf", &rates[n]) != 1) {
	fprintf(stderr, "Error reading <processor.rates>\n");
	signal_failure(1);
      }
      norm += rates[n];
    }
    fclose(in);

    for (n=0; n<numprocs; n++) {

      rates[n] /= norm;

      if (n == 0)
	nbodies_table[n] = nbodies_index[n] = round_up(rates[n] * nbodies_tot);
      else {
	if (n < numprocs-1)
	  nbodies_index[n] = round_up(rates[n] * nbodies_tot) + 
	    nbodies_index[n-1];
	else
	  nbodies_index[n] = nbodies_tot;
      
	nbodies_table[n] = nbodies_index[n] - nbodies_index[n-1];
      }

    }

    out = fopen("current.processor.rates", "w");
    if (out != NULL) {
      fprintf(out, "# %10s   %10s   %10s   %10s   %10s\n",
	      "Norm rate", "Raw rate", "Delta rate", "Index", "Current #");
      fprintf(out, "# %10s   %10s   %10s   %10s   %10s\n",
		"---------", "--------", "----------", "--------", "---------");
      for (n=0; n<numprocs; n++)
	fprintf(out, "  %10.7f   %10.3f   %10.7f   %10d   %10d\n",
		rates[n], 
		rates[n]*norm, 
		1.0 - rates[n]*nbodies_tot/nbodies_table[n],
		nbodies_index[n],
		nbodies_table[n]);
      fclose(out);
    }

  }

  MPI_Allreduce(&fail, &fail0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (fail0) {
    MPI_Finalize();
    exit(-1);
  }

  MPI_Bcast(nbodies_index, numprocs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nbodies_table, numprocs, MPI_INT, 0, MPI_COMM_WORLD);

}

void recurse_norm(double maxrate, double* r);


void normalize_rates(void)
{
  int n;
  int insane=0;
  double norm1;
				/* Particle fraction can't be more     */
				/* than this without overwriting array */

  maxrate = (double)nbodmax/(double)nbodies;

				/* Copy and normalize raw rates */
  norm1 = 0.0;
  for (n=0; n<numprocs; n++) {
    trates[n] = orates[n];
				/* Sanity check */
    if (trates[n]<0.0) {
      trates[n] = 0.0;
      insane=1;
    }

    norm1 += trates[n];
  }
  
  for (n=0; n<numprocs; n++)
    trates[n] /= norm1;

				/* Fix any overflows */
  recurse_norm(maxrate, trates);

				/* Test */
  norm1 = 0.0;
  for (n=0; n<numprocs; n++)
    norm1 += trates[n];
  
  if (fabs(norm1-1.0) > 1.0e-10) {
    fprintf(stderr, "\nNormalize_rates: norm out of bounds [%f]\n", norm1);
    for (n=0; n<numprocs; n++) {
      fprintf(stderr, "%3d:  %f  %f  %d\n", 
	      n, trates[n], orates[n], nbodies_table[n]);
      trates[n] /= norm1;
    }
  }

  if (insane) {
    fprintf(stderr, "\nNormalize_rates: insane values\n");
    for (n=0; n<numprocs; n++)
      fprintf(stderr, "%3d:  %f  %f  %d\n", 
	      n, trates[n], orates[n], nbodies_table[n]);
  }

}

/*

  Crude load balancing routine

*/
void recompute_processor_rates(void)
{
  double newrate, receive, norm, max_duration, min_duration, duration;
  int n, redistribute=0;
  MPI_Status status;
  double begin_time;

				/* Get rates from all processors */
  if (myid==0) {

    begin_time = MPI_Wtime();

    norm = 0.0;

    orates[0] = (double)nbodies*nbalance/MPL_read_timer(1);
    
    for (n=1; n<numprocs; n++) {
      MPI_Recv(&receive, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD,
	       &status);

      orates[status.MPI_SOURCE] = receive;
    }

    normalize_rates();

				/* Get maximum rate */
    min_duration = 1.0e32;
    max_duration = 0.0;
    for (n=0; n<numprocs; n++) {
      duration = nbodies_table[n]/trates[n];

      min_duration = duration < min_duration ? duration : min_duration;
      max_duration = duration > max_duration ? duration : max_duration;
    }

				/* Compute rate tolerance */

    if (max_duration - min_duration > tbalance*min_duration) redistribute = 1;

  }
  else {
    
    newrate = (double)nbodies*nbalance/MPL_read_timer(1);

    MPI_Send(&newrate, 1, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);
  }

  MPI_Bcast ( &redistribute, 1, MPI_INT, 0, MPI_COMM_WORLD);
  

				/* Force redistribution for debugging */
#ifdef DEBUG
  redistribute = 1;
#endif

  if (redistribute) {
    
    if (myid == 0)
      for (n=0; n<numprocs; n++) rates[n] = trates[n];
    
    redistribute_particles();

  }

				/* Diagnostic output */
  if (myid == 0) {

    FILE* out;

      out = fopen("current.processor.rates", "w");
      if (out != NULL) {
	fprintf(out, "# %10s   %10s   %10s   %10s   %10s\n",
		"Norm rate", "Raw rate", "Delta rate", "Index", "Current #");
	fprintf(out, "# %10s   %10s   %10s   %10s   %10s\n",
		"---------", "--------", "----------", "--------", "---------");
	for (n=0; n<numprocs; n++)
	  fprintf(out, "  %10.7f   %10.3f   %10.7f   %10d   %10d\n",
		  trates[n], 
		  orates[n], 
		  1.0 - trates[n]*nbodies_tot/nbodies_table[n],
		  nbodies_index[n],
		  nbodies_table[n]);
	fclose(out);
      }

  }

}


void parallel_gather_coefficients(void)
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

void parallel_distribute_coefficients(void)
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


void compute_mfp(void)
{
  int i;

  for (i=1; i<=nbodies; i++) {
    mfp[i] = 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]) +
      mass[i]*(pot[i] + potext[i]);
  }
}


void test_mpi(void)
{
  int n;
  
  fprintf(stderr, "Process %d:    x=%f   y=%f   z=%f\n", myid,
	  x[nbodies], y[nbodies], z[nbodies]);
  fprintf(stderr, "Process %d:   vx=%f  vy=%f  vz=%f\n", myid,
	  vx[nbodies], vy[nbodies], vz[nbodies]);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(-10);

}


static double MPL_accum_time=0.0;
static double MPL_last_time=0.0;

void MPL_reset_timer(void)
{
  MPL_accum_time=0.0;
}

void MPL_start_timer(void)
{
  MPL_last_time = MPI_Wtime();
}

void MPL_stop_timer(void)
{
  double curtime = MPI_Wtime();

  MPL_accum_time += curtime - MPL_last_time;
  MPL_last_time = curtime;
}

double MPL_read_timer(int reset)
{
  double save = MPL_accum_time;
  if (reset) MPL_accum_time = 0.0;
  return save;
}

