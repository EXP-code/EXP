/*
  Output information about the system state to
  the log and body data files.
*/

#include "expand.h"

void out_diag(int);
void synchronize_velocity(int);

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void out_put(int n)
{

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid == 0) fprintf(stderr,"%d.",n); fflush(stderr);

  if (orbtrace && myid == norbid) orb_trace();

  if (n==0) {

    out_log(n);
    if (olist) out_list(n);
    if (tipsy) out_tipsy(n);
    if (diag && myid==0) out_diag(n);
    if (coef && myid==0) out_coef(n);
    if (cylinder && coefcyl && myid==0) out_coef_cyl(n);
    if (relx) out_relx(n);
    if (pcaout && myid==0) out_pca(n);
  }
  else if ( !(n%nlog) || !(n%nlist) || !(n%ntipsy) || !(n%ndiag) || !(n%ncoef) || !(n%nrelx) || !(n%npcaout) || !(n%nbalance) || !(n%nchkpt) || finish) {

    int not_sync=1;
    /*
    double begin_time = MPI_Wtime();
    */

    if (!(n%nlog)  || finish) {
      if (not_sync) {synchronize_velocity(1); not_sync = 0;}
      out_log(n);
    }

    if (olist) {
      if (!(n%nlist) || finish) {
	if (not_sync) {synchronize_velocity(1); not_sync = 0;}
	out_list(n);
      }
    }

    if (nchkpt) {
      if (!(n%nchkpt)) {
	if (not_sync) {synchronize_velocity(1); not_sync = 0;}
	out_chkpt();
      }
    }

    if (tipsy) {
      if (!(n%ntipsy) || finish) {
	if (not_sync) {synchronize_velocity(1); not_sync = 0;}
	out_tipsy(n);
      }
    }

    if (diag) {
      if (!(n%ndiag) || finish) {
	if (not_sync) {synchronize_velocity(1); not_sync = 0;}
	if (myid==0) out_diag(n);
      }
    }

    if (coef) {
      if (!(n%ncoef) || finish) {
	if (myid==0) out_coef(n);
      }
    }
    
    if (cylinder && coefcyl) {
      if (!(n%ncoefcyl) || finish) {
	if (myid==0) out_coef_cyl(n);
      }
    }
    
    if (relx) {
      if (!(n%nrelx) || finish) {
	if (not_sync) {synchronize_velocity(1); not_sync = 0;}
	out_relx(n);
      }
    }

    if (pcaout) {
      if (!(n%npcaout) || finish) {
	if (myid==0) out_pca(n);
      }
    }

    if (!not_sync) synchronize_velocity(0);

    if (balance) {
      if (!(n%nbalance)) {
	recompute_processor_rates();
      }
    }

  }

}




