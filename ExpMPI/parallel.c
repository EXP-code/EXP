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

  if (myid == 0) {

    if ( (in = fopen("processor.rates", "r")) == NULL ) {
      fprintf(stderr, "setup: Error opening <processor.rates> . . . quit\n");
      signal_failure(1);
    }

    rates = (double *)malloc( slaves * sizeof(double) );
    if (rates==NULL) {
      fprintf(stderr, "setup: Error inializing rates\n");
      signal_failure(1);
    }
    rates--;

    orates = (double *)malloc( slaves * sizeof(double) );
    if (orates==NULL) {
      fprintf(stderr, "setup: Error inializing orates\n");
      signal_failure(1);
    }
    orates--;

    trates = (double *)malloc( slaves * sizeof(double) );
    if (trates==NULL) {
      fprintf(stderr, "setup: Error inializing trates\n");
      signal_failure(1);
    }
    trates--;

    nbodies_index = (int *)malloc( numprocs * sizeof(int) );
    if (nbodies_index==NULL) {
      fprintf(stderr, "setup: Error inializing index\n");
      signal_failure(1);
    }

    nbodies_table = (int *)malloc( slaves * sizeof(int) );
    if (nbodies_table==NULL) {
      fprintf(stderr, "setup: Error inializing table\n");
      signal_failure(1);
    }
    nbodies_table--;


    for (n=1; n<=slaves; n++) {
      if (fscanf(in, "%lf", &rates[n]) != 1) {
	fprintf(stderr, "Error reading <processor.rates>\n");
	signal_failure(1);
      }
      norm += rates[n];
    }
    fclose(in);

    nbodies_index[0] = 0;
    for (n=1; n<=slaves; n++) {

      rates[n] /= norm;

      if (n < slaves)
	nbodies_index[n] = round_up(rates[n] * nbodies) + nbodies_index[n-1];
      else
	nbodies_index[n] = nbodies;
      
      nbodies_table[n] = nbodies_index[n] - nbodies_index[n-1];

    }

    out = fopen("current.processor.rates", "w");
    if (out != NULL) {
      for (n=1; n<=slaves; n++)
	fprintf(out, "%f  %f  %f  %d\n", 
		rates[n], rates[n]*norm, 1.0 - rates[n]*nbodies/nbodies_table[n],
		nbodies_table[n]);
      fclose(out);
    }

    signal_failure(0);		/* OK! */

  } 
  else {

    int fail;			/* Will return zero if rate setup is ok */

    MPI_Bcast(&fail, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (fail) {
      MPI_Finalize();
      exit(-1);
    }
  }

}

				/* NR routine */
void indexx(int n, double *arrin, int *indx);

				/* Recursive normalization of         */
				/* rates to optimze n-body throughput */

void recurse_norm(double excess, double* r, int nleft)
{
				/* Waste an entry to avoid offsets */
  int *indx = (int *)malloc((nleft+1)*sizeof(int));
  int n;
  double norm2;

  if (indx == 0) {
    fprintf(stderr, "recurse_norm: error allocating temp storage\n");
    exit(-1);
  }

				/* Distribute excess among remainder */

  norm2 = 0.0;
  for (n=1; n<=nleft; n++) norm2 += r[n];
  for (n=1; n<=nleft; n++) r[n] += excess * r[n]/norm2;


				/* Get rank order */
  indexx(nleft, r, indx);

				/* Check for rate overflow */
    
  excess = 0.0;

  for (n=nleft; n>0; n--) {

    if (r[indx[n]] >= maxrate) {
      excess += r[indx[n]] - maxrate;
      r[indx[n]] = maxrate;
    }
    else if (excess > 0.0) {

      if (n == 1)
	r[indx[n]] += excess;
      else {
	int nn;
	double * rem = (double *)malloc((n+1)*sizeof(double));

	if (rem == 0) {
	  fprintf(stderr, "recurse_norm: error allocating temp storage\n");
	  exit(-1);
	}

	for (nn=1; nn<=n; nn++) rem[nn] = r[indx[nn]];
	recurse_norm(excess, rem, n);
	for (nn=1; nn<=n; nn++) r[indx[nn]] = rem[nn];

	free(rem);
      }
      
      free(indx);
      return;
    }
  }

  free(indx);
}

void normalize_rates(void)
{
  int n;
  int insane=0;
  double norm1, excess;

				/* Particle fraction can't be more     */
				/* than this without overwriting array */

  maxrate = (double)nbodmax/(double)nbodies;

				/* Copy and normalize raw rates */
  norm1 = 0.0;
  for (n=1; n<=slaves; n++) {
    trates[n] = orates[n];
				/* Sanity check */
    if (trates[n]<0.0) {
      trates[n] = 0.0;
      insane=1;
    }

    norm1 += trates[n];
  }
  
  for (n=1; n<=slaves; n++)
    trates[n] /= norm1;

				/* Check for storage overfill */
  excess = 0.0;
  recurse_norm(0.0, trates, slaves);

				/* Test */
  norm1 = 0.0;
  for (n=1; n<=slaves; n++)
    norm1 += trates[n];
  
  if (fabs(norm1-1.0) > 1.0e-10) {
    fprintf(stderr, "\nNormalize_rates: norm out of bounds [%f]\n", norm1);
    for (n=1; n<=slaves; n++) {
      fprintf(stderr, "%3d:  %f  %f  %d\n", 
	      n, trates[n], orates[n], nbodies_table[n]);
      trates[n] /= norm1;
    }
  }

  if (insane) {
    fprintf(stderr, "\nNormalize_rates: insane values\n");
    for (n=1; n<=slaves; n++)
      fprintf(stderr, "%3d:  %f  %f  %d\n", 
	      n, trates[n], orates[n], nbodies_table[n]);
  }

}

void distribute_particles(void)
{
  int n;
  MPI_Status status0;
  MPI_Status status[12];
  MPI_Request * request = (MPI_Request*)malloc((slaves*12)*sizeof(MPI_Request));

  if (request == 0) {
    fprintf(stderr, "distribute_particles: error allocating temp storage\n");
    exit(-1);
  }


#ifdef MPE_PROFILE
  MPE_Log_event(1, myid, "send_p");
#endif

  /*=====================================*/
  /* Divide up particles into processors */
  /*=====================================*/

  if (myid == 0) {

    if (is_init) {
      MPI_Bcast(&restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&tnow, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    for (n=1; n<=slaves; n++) {
      MPI_Send(&nbodies_table[n], 1, MPI_INT, n, 10, MPI_COMM_WORLD);
      fprintf(stderr, "Process 0 sending Process %d: %d bodies\n", n, 
	      nbodies_table[n]);

      MPI_Isend(&component[nbodies_index[n]-nbodies_table[n]+1], 
		nbodies_table[n], MPI_INT, n, 18, 
		MPI_COMM_WORLD, &request[(n-1)*12 + 0]);

      MPI_Isend(&initial_mass[nbodies_index[n]-nbodies_table[n]+1], 
		nbodies_table[n], MPI_DOUBLE, n, 19, 
		MPI_COMM_WORLD, &request[(n-1)*12 + 1]);
      MPI_Isend(&mass[nbodies_index[n]-nbodies_table[n]+1], 
		nbodies_table[n], MPI_DOUBLE, n, 20, 
		MPI_COMM_WORLD, &request[(n-1)*12 + 2]);

      MPI_Isend(&x[nbodies_index[n]-nbodies_table[n]+1], 
		nbodies_table[n], MPI_DOUBLE, n, 21, 
		MPI_COMM_WORLD, &request[(n-1)*12 + 3]);
      MPI_Isend(&y[nbodies_index[n]-nbodies_table[n]+1],
		nbodies_table[n], MPI_DOUBLE, n, 22,
		MPI_COMM_WORLD, &request[(n-1)*12 + 4]);
      MPI_Isend(&z[nbodies_index[n]-nbodies_table[n]+1],
		nbodies_table[n], MPI_DOUBLE, n, 23,
		MPI_COMM_WORLD, &request[(n-1)*12 + 5]);

      MPI_Isend(&vx[nbodies_index[n]-nbodies_table[n]+1],
		nbodies_table[n], MPI_DOUBLE, n, 24,
		MPI_COMM_WORLD, &request[(n-1)*12 + 6]);
      MPI_Isend(&vy[nbodies_index[n]-nbodies_table[n]+1], 
		nbodies_table[n], MPI_DOUBLE, n, 25,
		MPI_COMM_WORLD, &request[(n-1)*12 + 7]);
      MPI_Isend(&vz[nbodies_index[n]-nbodies_table[n]+1],
		nbodies_table[n], MPI_DOUBLE, n, 26,
		MPI_COMM_WORLD, &request[(n-1)*12 + 8]);

      MPI_Isend(&ax[nbodies_index[n]-nbodies_table[n]+1],
		nbodies_table[n], MPI_DOUBLE, n, 27,
		MPI_COMM_WORLD, &request[(n-1)*12 + 9]);
      MPI_Isend(&ay[nbodies_index[n]-nbodies_table[n]+1],
		nbodies_table[n], MPI_DOUBLE, n, 28,
		MPI_COMM_WORLD, &request[(n-1)*12 + 10]);
      MPI_Isend(&az[nbodies_index[n]-nbodies_table[n]+1],
		nbodies_table[n], MPI_DOUBLE, n, 29,
		MPI_COMM_WORLD, &request[(n-1)*12 + 11]);
    }


  }
  else {

    if (is_init) {
      MPI_Bcast(&restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
      fprintf(stderr, "Process %d:  restart=%d\n", myid, restart);

      MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      fprintf(stderr, "Process %d:  rmax=%f\n", myid, rmax);

      MPI_Bcast(&tnow, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      fprintf(stderr, "Process %d:  tnow=%f\n", myid, tnow);

      tpos = tvel = tnow;
    }

    MPI_Recv(&nbodies, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,
	     &status0);
    fprintf(stderr, "Process %d: nbodies=%d\n", myid, nbodies);
    
    MPI_Irecv(&component[1], nbodies, MPI_INT, MPI_ANY_SOURCE, 18,
	      MPI_COMM_WORLD, &request[0]);

    MPI_Irecv(&initial_mass[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 19,
	      MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(&mass[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 20,
	      MPI_COMM_WORLD, &request[2]);

    MPI_Irecv(&x[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 21,
	      MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(&y[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 22,
	      MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(&z[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 23,
	      MPI_COMM_WORLD, &request[5]);

    MPI_Irecv(&vx[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 24,
	      MPI_COMM_WORLD, &request[6]);
    MPI_Irecv(&vy[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 25,
	      MPI_COMM_WORLD, &request[7]);
    MPI_Irecv(&vz[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 26,
	      MPI_COMM_WORLD, &request[8]);

    MPI_Irecv(&ax[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 27,
	      MPI_COMM_WORLD, &request[9]);
    MPI_Irecv(&ay[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 28,
	      MPI_COMM_WORLD, &request[10]);
    MPI_Irecv(&az[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 29,
	      MPI_COMM_WORLD, &request[11]);

    MPI_Waitall(12, request, status);

  }

#ifdef MPE_PROFILE
  MPE_Log_event(2, myid, "recv_p");
#endif

  if (scatter) distribute_mfp();
  if (!is_init && relx) distribute_relx();

				/* Look for point mass particles */
  make_pointmass_list();

  free(request);
}

/*

  Crude load balancing routine

*/
void recompute_processor_rates(int *not_gathered)
{
  double newrate, receive, norm, max_duration, min_duration, duration;
  int n, redistribute=0;
  MPI_Status status;
  double begin_time;

				/* Get rates from all processors */
  if (myid==0) {

    begin_time = MPI_Wtime();

    norm = 0.0;

    for (n=1; n<=slaves; n++) {
      MPI_Recv(&receive, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD,
	       &status);

      orates[status.MPI_SOURCE] = receive;
    }

    normalize_rates();

				/* Get maximum rate */
    min_duration = 1.0e32;
    max_duration = 0.0;
    for (n=1; n<=slaves; n++) {
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
  

  if (redistribute) {
    
    if (*not_gathered) {
      gather_particles();
      *not_gathered = 0;
    }

    if (myid == 0) {

      for (n=1; n<=slaves; n++) rates[n] = trates[n];

				/* Recompute indices and particle numbers */
      nbodies_index[0] = 0;
      for (n=1; n<=slaves; n++) {

	if (n < slaves)
	  nbodies_index[n] = round_up(rates[n] * nbodies) + nbodies_index[n-1];
	else
	  nbodies_index[n] = nbodies;
	
	nbodies_table[n] = nbodies_index[n] - nbodies_index[n-1];
	
      }

    }

    distribute_particles();

  }

				/* Diagnostic output */
  if (myid == 0) {

    FILE* out;

      out = fopen("current.processor.rates", "w");
      if (out != NULL) {
	for (n=1; n<=slaves; n++)
	  fprintf(out, "%f  %f  %f  %d\n", 
		  trates[n], orates[n], 1.0 - trates[n]*nbodies/nbodies_table[n],
		  nbodies_table[n]);
	fclose(out);
      }

  }

}

void gather_particles(void)
{
  int n;
  MPI_Status status;


#ifdef MPE_PROFILE
  MPE_Log_event(3, myid, "b_gather_p");
#endif

  /*===============================*/
  /* Collect particles from slaves */
  /*===============================*/

  if (myid == 0) {

    for (n=1; n<=slaves; n++) {

				/* Tell slave to start sending */
      MPI_Send(&n, 1, MPI_INT, n, 9, MPI_COMM_WORLD);


      MPI_Recv(&component[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_INT, n, 18, MPI_COMM_WORLD, &status);

      MPI_Recv(&mass[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 20, MPI_COMM_WORLD, &status);

      MPI_Recv(&x[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 21, MPI_COMM_WORLD, &status);

      MPI_Recv(&y[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 22, MPI_COMM_WORLD, &status);

      MPI_Recv(&z[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 23, MPI_COMM_WORLD, &status);

      MPI_Recv(&vx[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 24, MPI_COMM_WORLD, &status);

      MPI_Recv(&vy[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 25, MPI_COMM_WORLD, &status);

      MPI_Recv(&vz[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 26, MPI_COMM_WORLD, &status);

      MPI_Recv(&ax[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 27, MPI_COMM_WORLD, &status);
	 
      MPI_Recv(&ay[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 28, MPI_COMM_WORLD, &status);

      MPI_Recv(&az[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 29, MPI_COMM_WORLD, &status);
      
      MPI_Recv(&pot[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 30, MPI_COMM_WORLD, &status);

      MPI_Recv(&potext[nbodies_index[n]-nbodies_table[n]+1],
	       nbodies_table[n], MPI_DOUBLE, n, 31, MPI_COMM_WORLD, &status);

    }

  }
  else {

    int nn;

				/* Wait for message from master before sending */

    MPI_Recv(&nn, 1, MPI_INT, 0, 9, MPI_COMM_WORLD, &status);

				/* Ok, go! */

    MPI_Send(&component[1], nbodies, MPI_INT, 0, 18, MPI_COMM_WORLD);

    MPI_Send(&mass[1], nbodies, MPI_DOUBLE, 0, 20, MPI_COMM_WORLD);

    MPI_Send(&x[1], nbodies, MPI_DOUBLE, 0, 21, MPI_COMM_WORLD);

    MPI_Send(&y[1], nbodies, MPI_DOUBLE, 0, 22, MPI_COMM_WORLD);

    MPI_Send(&z[1], nbodies, MPI_DOUBLE, 0, 23, MPI_COMM_WORLD);

    MPI_Send(&vx[1], nbodies, MPI_DOUBLE, 0, 24, MPI_COMM_WORLD);

    MPI_Send(&vy[1], nbodies, MPI_DOUBLE, 0, 25, MPI_COMM_WORLD);

    MPI_Send(&vz[1], nbodies, MPI_DOUBLE, 0, 26, MPI_COMM_WORLD);

    MPI_Send(&ax[1], nbodies, MPI_DOUBLE, 0, 27, MPI_COMM_WORLD);

    MPI_Send(&ay[1], nbodies, MPI_DOUBLE, 0, 28, MPI_COMM_WORLD);

    MPI_Send(&az[1], nbodies, MPI_DOUBLE, 0, 29, MPI_COMM_WORLD);

    MPI_Send(&pot[1], nbodies, MPI_DOUBLE, 0, 30, MPI_COMM_WORLD);

    MPI_Send(&potext[1], nbodies, MPI_DOUBLE, 0, 31, MPI_COMM_WORLD);
  }

#ifdef MPE_PROFILE
  MPE_Log_event(4, myid, "e_gather_p");
#endif
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

/* Send energies to slaves for keeping track of relaxation rate */

void distribute_relx(void)
{

  int i, n, nn;
  static int firstime=1;
  MPI_Status status0;
  MPI_Status *status = (MPI_Status*)malloc((slaves)*sizeof(MPI_Status));
  MPI_Request *request = (MPI_Request*)malloc((slaves)*sizeof(MPI_Request));

  if (request == 0 || slaves == 0) {
    fprintf(stderr, "distribute_relx: error allocating temp storage\n");
    exit(-1);
  }

#ifdef MPE_PROFILE
  MPE_Log_event(15, myid, "b_relx");
#endif

  if (firstime) {

    if (myid==0)
      esave = (double *)malloc((unsigned) nbodies*sizeof(double));
    else
      esave = (double *)malloc((unsigned) nbodmax*sizeof(double));

    if (!esave) {
      fprintf(stderr,"Processor %d: couldn't allocate energy storage vector\n",
	      myid);
      exit(-1);
    }
    esave -= 1;

				/* Get fiducial energy from slaves */
    if (myid == 0) {

      for (n=1; n<=slaves; n++) {
				/* Tell slave to start sending */
	MPI_Send(&n, 1, MPI_INT, n, 9, MPI_COMM_WORLD);

	MPI_Recv(&esave[nbodies_index[n]-nbodies_table[n]+1],
		 nbodies_table[n], MPI_DOUBLE, n, 32, MPI_COMM_WORLD, &status0);

      }
    }
    else {

      for (i=1; i<=nbodies; i++) {
	esave[i] = 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]) +
	  mass[i]*(pot[i] + potext[i]);
      }

				/* Wait for message from master before sending */

      MPI_Recv(&nn, 1, MPI_INT, 0, 9, MPI_COMM_WORLD, &status0);

				/* Ok, go! */

      MPI_Send(&esave[1], nbodies, MPI_DOUBLE, 0, 32, MPI_COMM_WORLD);

    }

    
    firstime = 0;

  }

      
				/* Redistribute */

  if (myid == 0) {

    for (n=1; n<=slaves; n++) {
      
      MPI_Isend(&esave[nbodies_index[n]-nbodies_table[n]+1], 
		nbodies_table[n], MPI_DOUBLE, n, 32, MPI_COMM_WORLD,
		&request[n-1]);
    }

  }
  else {

    MPI_Irecv(&esave[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 32,
	     MPI_COMM_WORLD, &request[0]);
  }


#ifdef MPE_PROFILE
  MPE_Log_event(16, myid, "e_relx");
#endif

  free(status);
  free(request);

}


/* Send int(rho*v) to slaves for keeping track of relaxation rate */

void distribute_mfp(void)
{

  int i, n, nn;
  static int firstime=1;
  MPI_Status status0;
  MPI_Status *status = (MPI_Status*)malloc((slaves)*sizeof(MPI_Status));
  MPI_Request *request = (MPI_Request*)malloc((slaves)*sizeof(MPI_Request));

  if (request == 0 || slaves == 0) {
    fprintf(stderr, "distribute_mfp: error allocating temp storage\n");
    exit(-1);
  }

  if (firstime) {

    if (myid==0)
      mfp = (double *)malloc((unsigned) nbodies*sizeof(double));
    else
      mfp = (double *)malloc((unsigned) nbodmax*sizeof(double));

    if (!mfp) {
      fprintf(stderr,"Processor %d: couldn't allocate mfp storage vector\n",
	      myid);
      exit(-1);
    }
    mfp -= 1;
    
				/* Get fiducial energy from slaves */
    if (myid == 0) {

      for (n=1; n<=slaves; n++) {
				/* Tell slave to start sending */
	MPI_Send(&n, 1, MPI_INT, n, 9, MPI_COMM_WORLD);

	MPI_Recv(&mfp[nbodies_index[n]-nbodies_table[n]+1],
		 nbodies_table[n], MPI_DOUBLE, n, 33, MPI_COMM_WORLD, &status0);

      }
    }
    else {

      for (i=1; i<=nbodies; i++) {
	mfp[i] = 0.5*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]) +
	  mass[i]*(pot[i] + potext[i]);
      }

				/* Wait for message from master before sending */

      MPI_Recv(&nn, 1, MPI_INT, 0, 9, MPI_COMM_WORLD, &status0);

				/* Ok, go! */
      
      MPI_Send(&mfp[1], nbodies, MPI_DOUBLE, 0, 33, MPI_COMM_WORLD);

    }

    
    firstime = 0;

  }

      
				/* Redistribute */

  if (myid == 0) {

    for (n=1; n<=slaves; n++) {
      
      MPI_Isend(&mfp[nbodies_index[n]-nbodies_table[n]+1], 
		nbodies_table[n], MPI_DOUBLE, n, 33, MPI_COMM_WORLD,
		&request[n-1]);
    }

  }
  else {

    MPI_Irecv(&mfp[1], nbodies, MPI_DOUBLE, MPI_ANY_SOURCE, 33,
	      MPI_COMM_WORLD, &request[0]);
  }

  free(status);
  free(request);
}


void test_mpi(void)
{
  int n;
  
  if (myid > 0) {
    fprintf(stderr, "Process %d:    x=%f   y=%f   z=%f\n", myid,
	    x[nbodies], y[nbodies], z[nbodies]);
    fprintf(stderr, "Process %d:   vx=%f  vy=%f  vz=%f\n", myid,
	    vx[nbodies], vy[nbodies], vz[nbodies]);
  }
  else {
    for (n=1; n<=slaves; n++) {
      fprintf(stderr, "Process %d/0:  x=%f   y=%f   z=%f\n", n,
	     x[nbodies_index[n]], y[nbodies_index[n]], z[nbodies_index[n]]);
      fprintf(stderr, "Process %d/0: vx=%f  vy=%f  vz=%f\n", n,
	     vx[nbodies_index[n]], vy[nbodies_index[n]], vz[nbodies_index[n]]);
    }
    
  }

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

