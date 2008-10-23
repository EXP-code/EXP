/*
  Get the whole thing going by:
  -----------------------------
  1) reading in initial phase space
  2) computing potential initially
  3) dumping first phase-space and log file
  4) write parameter file
*/

#include <expand.h>
#include <ExternalCollection.H>
#include <OutputContainer.H>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void begin_run(void)
{

  //===================================
  // Make the Vector class mutex
  //===================================
  
  int errcode = pthread_mutex_init(&mem_lock, NULL);
  if (errcode) {
    cerr << "Process " << myid 
	 << ": failed to make the Vector class memory lock, errcode=" 
	 << errcode << endl;
    MPI_Abort(MPI_COMM_WORLD, 115);
  }

  //===================================
  // Make the kick/drift thread vectors
  //===================================
  
  posvel_data = vector<thrd_pass_posvel>(nthrds);
  posvel_thrd = vector<pthread_t>       (nthrds);

  //==============================
  // Initialize multistepping
  //==============================

  initialize_multistep();

  //===================================
  // Initialize phase-space components 
  //===================================

  comp.initialize();

  //===================================
  // Initialize external forces, if any
  //===================================

  external.initialize();
  
  //===============================
  // Compute initial acceleration
  //===============================

  initializing = true;

  if (multistep) {
    sync_eval_multistep();	// Use last coefficient evaluation
    for (int M=0; M<=multistep; M++) {
      comp.compute_expansion(M);
    }
  } else {
    comp.multistep_reset();
  }

  comp.compute_potential(0);
  //                     ^
  //                     |
  //   All time levels---/

  //==============================
  // Initialize multistep levels
  //==============================
  adjust_multistep_level(true);
  //                     ^
  //                     |
  // Do all particles---/

  // Then recompute coefficients . . . 
  //
  if (multistep) {
    sync_eval_multistep();	// Use last coefficient evaluation
    for (int M=0; M<=multistep; M++) {
      comp.compute_expansion(M);
    }
  } else {
    comp.multistep_reset();
  }

  comp.compute_potential(0);

  initializing = false;

  //===================================
  // Initialize output routines
  //===================================

  output.initialize();
  output.Run(0);

  //======================
  // Write parameter file 
  //======================

  write_parm();

}
